<?php
$sMiSeek = "/beegfs/group_dv/software/source/miRNASearch/mirseek/Release/mirseek";
$sLibDef = "libdef.txt";
$sMiProbes = "mrna_probes.fa";
$nMaxDiff = 0;
$nMaxMatchPerRead = 100;

if (count($argv) >=4) {
	$nMaxDiff =intval($argv[3]) ;
}

if (count($argv) >= 5) {
	$sLibDef = $argv[4];
}

$sDumpDir = "out_maxdiff$nMaxDiff";
$nReadLimit = 1e9; 

$nTotalParts = intval($argv[1]);
$nThisPart = intval($argv[2]);

$arrSamples = fnLoadSamples($sLibDef);
$arrProbes = fnLoadProbes($sMiProbes);

//print_r($arrProbes);
exec("mkdir -p $sDumpDir");
$nJob = -1;
foreach($arrSamples as $sTaxon => $arrReads) {


	foreach($arrProbes as $sProbeID => $sProbeSeq) {
		if ( (++$nJob) %  $nTotalParts != $nThisPart ) {
			continue;
		}

		echo("Job $nJob: $sTaxon $sProbeSeq ...\n");

		$sDumpTaxonDir = "$sDumpDir/$sTaxon";

		if (fnTestOutput("$sDumpTaxonDir/count_for_$sProbeSeq.txt")) {
			echo("$sTaxon $sProbeSeq done, skip");
			continue;
		}

		exec("mkdir -p $sDumpTaxonDir");

		$h1 = popen("zcat -f $arrReads[0]" , 'r');
		$h2 = popen("zcat -f $arrReads[1]" , 'r');

		$arrBuffer = array();

		$descriptorspec = array(
		   0 => array("pipe", "r"),  // stdin is a pipe that the child will read from
		   1 => array("pipe", "w"),  // stdout is a pipe that the child will write to
		   2 => array("pipe", "w") // stderr is a file to write to
		);

		$hMirSeek = proc_open("$sMiSeek", $descriptorspec, $arrPipes);
		$hDump1 = popen("gzip -c > $sDumpTaxonDir/hit_reads_for_$sProbeSeq.R1.fq.gz", 'w');
		$hDump2 = popen("gzip -c > $sDumpTaxonDir/hit_reads_for_$sProbeSeq.R2.fq.gz", 'w');
		$hO = fopen("$sDumpTaxonDir/count_for_$sProbeSeq.txt", 'w');

		fwrite($arrPipes[0], "$nMaxDiff\n");
		fwrite($arrPipes[0], "$nMaxMatchPerRead\n");
		fwrite($arrPipes[0], "$sProbeSeq\n");

		$nTotalHits = 0;
		$nTotalHit1 = 0;
		$nTotalHit2 = 0;
		$nTotalDump = 0;
		$nProccessedReads = 0;
		$nProccessedBases = 0;
		while(true) {
			if ($nProccessedReads > $nReadLimit ) {
				break;
			}
			$sLn1 = fgets($h1);
			$sLn2 = fgets($h2);

			if ($sLn1 === false || $sLn2 === false) break;

			if ($sLn1[0] != '@' || $sLn2[0] != '@') {
				die("First line unexpected, $arrReads[0], $arrReads[1]\n$sLn1\n$sLn2\n");
			}
			$arrBuffer = array(array('read_id'=>trim($sLn1) , 'seq' => '', 'qual' => ''), array('read_id'=>trim($sLn2) , 'seq' => '', 'qual' => ''));

			$sLn1 = fgets($h1);
			$sLn2 = fgets($h2);

			if ($sLn1 === false || $sLn2 === false) {
				echo("Warning: $arrReads[0], $arrReads[1] ended prematurely\n");
				break; 
			}
			$arrBuffer[0]['seq'] = trim($sLn1);
			$arrBuffer[1]['seq'] = trim($sLn2);

			$sLn1 = fgets($h1);
			$sLn2 = fgets($h2);

			if ($sLn1 === false || $sLn2 === false) {
				echo("Warning: $arrReads[0], $arrReads[1] ended prematurely\n");
				break; 
			}

			if ($sLn1[0] != '+' || $sLn2[0] != '+') {
				die("Third line unexpected, $arrReads[0], $arrReads[1]\n$sLn1\n$sLn2\n");
			}

			$sLn1 = fgets($h1);
			$sLn2 = fgets($h2);

			if ($sLn1 === false || $sLn2 === false) {
				echo("Warning: $arrReads[0], $arrReads[1] ended prematurely\n");
				break; 
			}
			$arrBuffer[0]['qual'] = trim($sLn1);
			$arrBuffer[1]['qual'] = trim($sLn2);

			fwrite($arrPipes[0], $arrBuffer[0]['seq']."\n");
			$nHit1 = fnReadNumFromStream($arrPipes[1]);
			fwrite($arrPipes[0], $arrBuffer[1]['seq']."\n");
			$nHit2 = fnReadNumFromStream($arrPipes[1]);

			if ( intval($nHit1 + $nHit2) >= 1) {
				$nTotalDump++;
				//echo("Hit1 $nHit1 Hit2 $nHit2 Write: R1 ".trim($arrBuffer[0]['read_id'])." ; R2 ".trim($arrBuffer[1]['read_id'])."\n");
				fwrite($hDump1 , $arrBuffer[0]['read_id']." READ1 $nHit1\n".$arrBuffer[0]['seq']."\n+\n".$arrBuffer[0]['qual']."\n");
				fwrite($hDump2 , $arrBuffer[1]['read_id']." READ2 $nHit2\n".$arrBuffer[1]['seq']."\n+\n".$arrBuffer[1]['qual']."\n");
			}

			$nTotalHits += $nHit1 + $nHit2;
			$nTotalHit1 += $nHit1;
			$nTotalHit2 += $nHit2;
			$nProccessedReads++;
			$nProccessedBases += strlen($arrBuffer[0]['seq']) + strlen($arrBuffer[1]['seq']);;
			if ($nProccessedReads % 1000000 ==0) {
				echo("Proccessed $nProccessedReads reads, $nProccessedBases bases, found $nTotalHits, written $nTotalDump probe-containing reads... \n");
			}

		}

		fwrite($arrPipes[0], "#STOP\n");
		$nProgramTotal = fnReadNumFromStream($arrPipes[1]);
		proc_close($hMirSeek);

		if ($nProgramTotal != $nTotalHits) {
			echo("Warning: total counts differ between host and php programs php $nTotalHits != host $nProgramTotal\n");
		} else {
			echo("Finished successfully\n");
		}

		fwrite($hO , "Taxon\tProbeID\tProbeSeq\tTotalHitR1\tTotalHitR2\tTotalHits\tProccessedReads\tProccessedBases\tWrittenReads\n");
		fwrite($hO , "$sTaxon\t$sProbeID\t$sProbeSeq\t$nTotalHit1\t$nTotalHit2\t$nTotalHits\t$nProccessedReads\t$nProccessedBases\t$nTotalDump\n");
		fclose($hO);
		pclose($hDump1);
		pclose($hDump2);
		echo("Total found read pairs: $nTotalDump\n");

		//die();
	}
	
}

function fnLoadSamples($s) {
	$h = fopen($s, 'r');
	$arrRet = array();
	while(false!==($sLn = fgets($h) )) {
		$sLn = trim($sLn);
		if ($sLn == '' || $sLn[0] == '#') continue;
		$arrF = explode("\t", $sLn);
		$sTaxon = $arrF[0];
		$sR1 = $arrF[3];
		$sR2 = $arrF[4];
		if (!file_exists($sR1)) {
			die("$sR1 not found\n");
		}

		if (!file_exists($sR2)) {
			die("$sR2 not found\n");
		}
		$arrRet[$sTaxon] = array($sR1, $sR2);
	}

	return $arrRet;
}

function fnLoadProbes($s) {
	$hIn = fopen($s, 'r');
	$sSeqName = "";
	$sSeq = "";
	$arrRet = array();
	while(true) {
		$sLn = fgets($hIn);
		        if (false === $sLn || (trim($sLn) != '' && $sLn[0] =='>')) {
		        if ($sSeq != '') {
		                $arrRet[str_replace(' ', '_', $sSeqName)] = $sSeq;
		        }
		        if (false === $sLn) {
		                break;
		        }
		}

		$sLn = trim($sLn);
		if ($sLn[0]=='>') {
		        $sSeqName = substr($sLn, 1);
		        $sSeq = "";
		        continue;
		}

		$sSeq .= strtoupper($sLn);
	}

	return $arrRet;
}

function fnReadNumFromStream($h) {
	$n = intval(trim(fgets($h)));
	/*while(false !== ($sLn=fgets($h)) ) {
		echo("extra line: $sLn");
	}*/

	return $n;
}

function fnTestOutput($sF) {
	$h = popen("cat $sF 2>/dev/null | wc -l", 'r');
	$nCount = trim(fgets($h));
	while(false!==fgets($h)) {};
	pclose($h);
	return $nCount=="2";

}

?>
