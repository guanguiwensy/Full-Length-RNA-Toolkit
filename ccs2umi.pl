#! /usr/bin/perl

open IN,"$ARGV[0]";
%mark1=();
while(<IN>){
	chomp;
	@a=();@a=split;
	@b=();@b=split(/\//,$a[0]);
	$tag="$b[0]$b[1]";
	$tag=~s/^>//;
	@c=();@c=split(/\;/,$a[1]);
	$c[0]=~s/strand=//;
	$c[2]=~s/polyAseen=//;
	$c[3]=~s/threeseen=//;
	$c[5]=~s/polyAend=//;
	$c[6]=~s/threeend=//;
	# print "$c[0]\t$c[1]\t$c[2]\t$c[3]\n";
	if($c[2]==1 and $c[3]==1){
		# print "$c[0]\t$c[2]\t$c[3]\t$c[5]\t$c[6]\n";
		$mark1{$tag}{strand}=$c[0];
		$mark1{$tag}{polyAend}=$c[5];
		$mark1{$tag}{threeend}=$c[6];
	}
}
close IN;

open IN,"$ARGV[1]";
%ccslen=();
while(<IN>){
	chomp;
	if($_=~/^>/){
		$name=$_;
		$name=~s/^>//;
		@a=();@a=split(/\//,$name);
		$tag="$a[0]$a[1]";
	}else{
		$len=length($_);
		$ccslen{$tag}=$len;
	}
}
close IN;

open IN,"$ARGV[2]";
%primer=();%seqstat=();
while(<IN>){
	chomp;
	if ($_=~/^m/){
		@a=();@a=split;
		@b=();@b=split(/\//,$a[0]);
		$tag="$b[0]$b[1]";
		if(defined $mark1{$tag}){
			if($mark1{$tag}{strand} eq "+"){
				$pos="$a[6]\_$a[7]";
				$len=$a[8]-$a[9]+1;
				if($a[8]>$a[9]){
					# print "$tag\t$ccslen{$tag}\t$mark1{$tag}{strand}\t$mark1{$tag}{threeend}\t$a[6]\t$a[7]\t$len\n";
					if(!defined $primer{$tag}){
						$dis=0;$dis=abs($mark1{$tag}{threeend}-$a[6]);
						$primer{$tag}[0]=$dis;
						$primer{$tag}[1]=$pos;
					}else{
						$dis=0;$dis=abs($mark1{$tag}{threeend}-$a[6]);
						if($dis<$primer{$tag}[0]){
							$primer{$tag}[0]=$dis;
							$primer{$tag}[1]=$pos;
						}
					}
					# $primer{$tag}[]=$pos if(!defined $primer{$tag});
				}
			}elsif($mark1{$tag}{strand} eq "-"){
				$pos="$a[6]\_$a[7]";
				$len=$a[9]-$a[8]+1;
				if($a[8]<$a[9]){
					# print "$tag\t$ccslen{$tag}\t$mark1{$tag}{strand}\t$mark1{$tag}{threeend}\t$a[6]\t$a[7]\t$len\n";
					if(!defined $primer{$tag}){
						$dis=0;$dis=abs($mark1{$tag}{threeend}-$a[7]);
						$primer{$tag}[0]=$dis;
						$primer{$tag}[1]=$pos;
					}else{
						$dis=0;$dis=abs($mark1{$tag}{threeend}-$a[7]);
						if($dis>$primer{$tag}[0]){
							$primer{$tag}[0]=$dis;
							$primer{$tag}[1]=$pos;
						}
					}
					# $primer{$tag}=$pos if(!defined $primer{$tag});
				}
			}
		}
	}
}
close IN;

open IN,"$ARGV[1]";
%ercc_umi=();
open OUT,">$ARGV[3]/ccs2umi.xls";
print OUT  "CCS_ID\tCCS_Len\tPrimer_Start\tPrimer_End\tUMI\tpolyA/T\n";
while (<IN>){
	chomp;
	if ($_=~/^>/){
		$name=$_;
		$name=~s/^>//;
		@a=();@a=split(/\//,$name);
		$tag="$a[0]$a[1]";
	}else{
		$seq="";$seq=$_;
		$mark1=0;$mark2=0;
		if ($mark1{$tag}{strand} eq "-"){
			@b=();@b=split(/\_/,$primer{$tag}[1]);
			$barcode=substr($seq,$b[1],6);
			$polyA=substr($seq,$b[1]+6,30);
			$reccomp_barcode=reverse($barcode);
			$reccomp_barcode=~tr/ACGTacgt/TGCAtgca/;
			print OUT "$name\t$ccslen{$tag}\t$b[0]\t$b[1]\t$reccomp_barcode\t$polyA\n" if ($barcode ne "");
			# print "$name\t$ccslen{$tag}\t$b[0]\t$b[1]\t$barcode\t$polyA\n" if ($barcode ne "");
		}
		if ($mark1{$tag}{strand} eq "+"){
			@b=();@b=split(/\_/,$primer{$tag}[1]);
			$barcode=substr($seq,$b[0]-7,6);
			$polyA=substr($seq,$b[0]-37,30);
			print OUT "$name\t$ccslen{$tag}\t$b[0]\t$b[1]\t$barcode\t$polyA\n" if ($barcode ne "");
		}
	}
}
close IN;

