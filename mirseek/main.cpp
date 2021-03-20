#include <stdio.h>
#include <cstring>
#include <ctype.h>
#include <string>
#include <iostream>
#include "main.h"
#include "edlib.h"
using namespace std;

int main(int argc, char **argv)
{
	//char query[] = "uAAGUGCUUCUCUUUGGGGUUG";
	//char target[] = "UAAGAGCUUCUCUUUGTGGUAGGCUCAGACCACCUCAAAGAAUCCACUGAUGUUUAAGUGCUUCUCUUUGGGGUUGUAAGUGCUUCUCUUUGGGGUUGUAAGUGCUUCUCUUUGGGGUUGUCUGCUCAUAAGUGCUUCUCUUUGGGGAAGUCUUGGUAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU";
	
	string sQuery;
	string sTarget;
	int nMaxMisMatches = 3;
	int nMaxHitPerRead = 100;
	cin >> nMaxMisMatches;
	cin >> nMaxHitPerRead;
	cin	>> sQuery;
	
	char * szQuery = new char[strlen(sQuery.c_str()) + 1];
	strcpy(szQuery, sQuery.c_str());
	strtoupper(szQuery);
	strU2T(szQuery); //convert U to T
		
	char * szQ_revcomp = new char[strlen(szQuery) + 1];
	fnRevComp(szQuery , szQ_revcomp);
		
	unsigned long long nTotal = 0;
	while(true) {
		cin >> sTarget;
		if (sTarget == "#STOP") {
			printf("%lld\n", nTotal);
			break;
		}
		int nHits = fnFindMatches(szQuery, sTarget.c_str(), nMaxMisMatches, nMaxHitPerRead);
		//printf("Totally %d hits found\n", nHits);
	
		if (strcmp(szQ_revcomp ,szQuery ) != 0 ) { // check if kmer is palindromic, if so no need to check rev comp.
			//printf("Rev comp: %s \n", szQ_revcomp);
			nHits += fnFindMatches(szQ_revcomp, sTarget.c_str(), nMaxMisMatches, nMaxHitPerRead);
			//printf("Totally %d hits found\n", nHits);		
		}
		printf("%d\n", nHits);
		nTotal += (unsigned long long)nHits;
	}
	delete[] szQuery;
	delete[] szQ_revcomp;

	return 0;
}

int fnFindMatches(const char * szQuery, const char * szTarget, int nMaxDiff, int nMaxMatches) {
	size_t nQLen = strlen(szQuery);
	size_t nTLen = strlen(szTarget);
	
	char * szQ = new char[nQLen+1];
	char * szT = new char[nTLen+1];
	strcpy(szQ, szQuery);
	strcpy(szT, szTarget);
	//strtoupper(szQ); //saves time
	strtoupper(szT);

	//int * arrEndPos = new int[nMaxMatches];
	int nMatches = 0;
	
	while(nMatches < nMaxMatches) {
		EdlibAlignResult result = edlibAlign(szQ, nQLen, szT, nTLen, 
									edlibNewAlignConfig(nMaxDiff, EDLIB_MODE_HW, EDLIB_TASK_LOC, additionalEqualities, nAdditionalEqualities) );
		if (result.status == EDLIB_STATUS_OK && result.numLocations>0) {
			//printf("Number of matches: %d\n", result.numLocations);
			int nLoc = 0;
			for(nLoc=0;nLoc<result.numLocations;nLoc++) {
				++nMatches;
				int nStart = result.startLocations[nLoc];
				int nEnd = result.endLocations[nLoc];
				//printf("at: %d - %d, diff = %d\n", nStart, nEnd, result.editDistance);
				//mask the found region
				memset(szT+nStart, '@', nEnd - nStart + 1);
				//printf("mask hit: %s\n", szT);
			}
		} else {
			edlibFreeAlignResult(result);
			break;
		}
		edlibFreeAlignResult(result);
	}
	
	delete[] szQ;
	delete[] szT;
	return (nMatches < nMaxMatches)? nMatches : nMaxMatches;
}

void strtoupper(char * s) {
  char * pOrig = s;
  while (*s) {
    *s = toupper((unsigned char) *s);
    s++;
  }
  s = pOrig;
}

void strU2T(char * s) {
  char * pOrig = s;
  while (*s) {
	  if(*s == 'U') {
		  *s = 'T';
	  }
    s++;
  }
  s = pOrig;
}

void fnRevComp(const char *szDNA, char * szOut) {
	int nLen = strlen(szDNA);
	int i = nLen - 1;
	char * pOrig = szOut;
	for(;i>=0;--i) {
		char nLetter = toupper(szDNA[i]);
		switch(nLetter) {
			case 'A': 
				*szOut = 'T'; break;
			case 'T': 
				*szOut = 'A'; break;
			case 'U': 
				*szOut = 'A'; break;
			case 'G': 
				*szOut = 'C'; break;
			case 'C': 
				*szOut = 'G'; break;
			case 'N': 
				*szOut = 'N'; break;
			default :
				*szOut = 'N'; break;
		}
		
		szOut++;
	}
	szOut = pOrig;
	szOut[nLen] = '\0';
}