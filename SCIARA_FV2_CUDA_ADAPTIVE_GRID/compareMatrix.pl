#! /usr/bin/perl
use File::Compare;

$ARGV[0] ||  die "path matrix 1 not found\n";
$ARGV[1] ||  die "path matrix 2 not found\n";
$VERBOSE=$ARGV[2];
if (compare($ARGV[0],$ARGV[1]) == 0) {
    print " $ARGV[0] and $ARGV[1] are equal\n";
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
}
	

