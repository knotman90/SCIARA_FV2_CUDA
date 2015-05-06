#! /usr/bin/perl
use File::Compare;
use File::Basename;


$ARGV[0] ||  die "*PATH matrix 1 not found*\n";
$ARGV[1] ||  die "*PATH matrix 2 not found*\n";


$PATH1 = $ARGV[0];
$PATH2 = $ARGV[1];

$VERBOSE=$ARGV[2];

print "Comparing:\n";
if($VERBOSE){

	print "	$ARGV[0]\n";
	print "	$ARGV[1]\n";
}else{
	my($file, $dir, $ext) = fileparse($PATH1);
	print basename "	$file\n";
	my($file, $dir, $ext) = fileparse($PATH2);
	print basename "	$file\n";
}


unless (-f $PATH1) { die "File 1 Doesn't Exist!\n"; }
unless (-f $PATH2) { die "File 2 Doesn't Exist!\n"; }


if (compare($PATH1,$PATH2) == 0) {
    die "!!!ATTENTION!!!\n $PATH1 and $PATH2 are THE SAME FILE\n";
}else
{
	
	open(FILE1,$ARGV[0]);
	open(FILE2,$ARGV[1]);
	$count=0;
	$sum = 0;
	$max = 0;
	$r=-1;
	$c=-1;
	while($line1=<FILE1>){
$c=0;
	$r=$r+1;
		$line2=<FILE2>;
		@valuesLine1 = split(/ /,$line1);
		@valuesLine2 = split(/ /,$line2);
		$size = @valuesLine1;
		for($i=0;$i<$size;$i++){
	$c=$c+1;
			if($valuesLine1[$i]!=$valuesLine2[$i]){
				$sum+=abs($valuesLine1[$i]-$valuesLine2[$i]);
				if(abs($valuesLine1[$i]-$valuesLine2[$i]) > 0 && $VERBOSE){
					print "($r,$c,$valuesLine1[$i],$valuesLine2[$i])\n";					
				}
				$count++;
				if(abs($valuesLine1[$i]-$valuesLine2[$i])>$max){
					$max = abs($valuesLine1[$i]-$valuesLine2[$i]);
				}
			}
		}
	}
	if($sum>0){
		$difference = $sum/$count;
		print "Num of differents cells: $count\n";
		print "Avarage error of differents cells: $difference\n";
		print "Max error: $max\n";
	}else{
		print "Matrices are equals\n";
	}
	close(FILE1);
	close(FILE2);

print "\n";
}
	

