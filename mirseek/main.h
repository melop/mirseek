#include "edlib.h"

//EdlibEqualityPair additionalEqualities[5] = {{'U', 'T'}, {'N', 'A'}, {'N', 'T'}, {'N', 'G'}, {'N', 'C'} };
//int nAdditionalEqualities = 5;

EdlibEqualityPair additionalEqualities[1] = {{'U', 'T'}};
int nAdditionalEqualities = 1;

int fnFindMatches(const char * szQuery, const char * szTarget, int nMaxDiff, int nMaxMatches);
void strtoupper(char * s) ;
void strU2T(char * s) ;
void fnRevComp(const char *szDNA, char * szOut) ;