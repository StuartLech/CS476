use strict;
use warnings;

my($data, @matrix, $i, $j, $x, $y, $val1, $val2);
my(@keys);
my($outfilename);

### BUNCH OF ADDED DECLARATIONS
my($fileinput,@temp);

#print "Enter the name of the mutation probability PAM matrix input file: ";
#$fileinput = <STDIN>;
#chomp $fileinput;
$fileinput = "PAM1-mutprob.txt";
if(!open(infile3, $fileinput)){
    print "error opening matrix file\n";
    exit;
}

# modified to ignore Header of substitution matrix files
$data = <infile3>; # ignore!

# read amino acids and build hash table
$data = <infile3>;
chomp $data;
@keys = split(/[,\s+]/, $data); # MODIFIED to remove , or spaces

# build hash table - map amino acids to array subscript
my %aminoacid;
for ($i=0; $i<20; $i++){
   $aminoacid{$keys[$i]} = $i;
}

my(@probmatrix);
my(@submatrix);
my(@tempmatrix);
my(@powxmatrix);

my(@rowsum);
my(@colsum);

#### HARDCODING AMINO ACID FREQUENCIES AFTER HAVING COMPUTED THEM FROM HIGH EXPONENTIATION PREVIOUSLY 333
my(@aminofreq);
@aminofreq = ( 0.08768333990119, 0.0405129182960021, 0.0408784631518651, 0.0477160345974603, 0.0324709539656211, 0.0378461268859548, 0.0504933695605074, 0.0898249006830963, 0.0328588505954496, 0.0357514442352249, 0.0852464099207531, 0.0791031344407513, 0.0148824394639692, 0.0410010190895639, 0.0515802694709073, 0.0697549720598532, 0.0583275704247605, 0.00931264523877659, 0.0317154088087089, 0.0630397292098708 );

# populate probability or substitution matrix
$j = 0;
while ($data = <infile3>){
   chomp $data;
   @temp = split(/[,\s+]/, $data);  # MODIFIED to remove , or spaces
   push @submatrix, [@temp];
   push @tempmatrix, [@temp];
   $j++;
}

my($probflag);
my($divergence);

print "Please enter the units of divergence to be computed:\n";
$divergence = int(<>);
#### Must raise the probability matrix to the power $divergence
#### This is the INEFFICIENT way to do so (hopefully to be commented out later)
my($powx,$zcount,$r,$c); # these are just count-like variables not used elsewhere (or overriding other declaration here otherwise)
my($cursum);

#compute row sum and column sum
#also compute total matrix sum
my($totsum);
$totsum = 0;
for ($r = 0; $r < 20; $r++) {
    $rowsum[$r] = 0; $colsum[$r] = 0;
    for ($c = 0; $c < 20; $c++) {
	$rowsum[$r] = $rowsum[$r] + $tempmatrix[$r][$c];
	$colsum[$r] = $colsum[$r] + $tempmatrix[$c][$r];
    }
    $totsum = $totsum + $rowsum[$r];
}

# Now get the actual STOCHASTIC matrix of probabilities:
for ($r = 0; $r < 20; $r++) {
    for ($c = 0; $c < 20; $c++) {
	$probmatrix[$r][$c] = $tempmatrix[$r][$c] / $colsum[$c];
    }
}

# now make the temp matrix also equal the mutation probabilities in the current form
for ($r = 0; $r < 20; $r++) {
    for ($c = 0; $c < 20; $c++) {
	$tempmatrix[$r][$c] = $probmatrix[$r][$c];
    }
}
    
## DEBUGGING!    
# Modified to Print out the probability matrix!
#    print "\n    ";
#   print "@keys\n";
#  for ($i = 0; $i < 20; $i++) {
#	print "$keys[$i]: ";
#	for ($j = 0; $j < 20; $j++) {
#	    print "$probmatrix[$i][$j] "; 
#	}
#	print "\n";
#   }

# matrix exponentiation (INEFFICIENT WAY!): now take the mutation probability matrix to the $divergence power
for ($powx = 2; $powx <= $divergence; $powx++) { # this is the outer loop iterating over the number of self-multiplications
    # remember: the tempmatrix is the same as the mutation probability matrix initially.  we will keep multiplying the current probmatrix by tempmatrix
    for ($r = 0; $r < 20; $r++) {
	for ($c = 0; $c < 20; $c++) {
	    $cursum = 0;
	    for ($zcount = 0; $zcount < 20; $zcount++) {
		$cursum = $cursum + $probmatrix[$r][$zcount]*$tempmatrix[$zcount][$c];
	    }
	    $powxmatrix[$r][$c] = $cursum;
	}
    }
    # now copy the $powxmatrix back into the $probmatrix
    for ($r = 0; $r < 20; $r++) {
	for ($c = 0; $c < 20; $c++) {
	    $probmatrix[$r][$c] = $powxmatrix[$r][$c];
	}
    }	
}

## FOR VERY HIGH DIVERGENCE, EVERY ENTRY OF A ROW i IS NOW THE FREQUENCY OF AMINO ACID i
print "\n    ";
print "FREQUENCIES: \n";
for ($i = 0; $i < 20; $i++) {
    print "$keys[$i]: ";
    for ($j = 0; $j < 20; $j++) {
	print "$probmatrix[$i][$j] "; 
    }
    print "\n";
}

$x = <STDIN>;
    
    
## NOW WE HAVE THE CORRECT MUTATION PROBABILITIES!
## TO OBTAIN THE ACTUAL SUBSTITUTION SCORES USE THE FOLLOWING EQUATION
## substitution score[x][y] = 10*log( Prob[x][y] / Freq[y])
for ($r = 0; $r < 20; $r++) {
    for ($c = 0; $c < 20; $c++) {
	#	    $submatrix[$r][$c] = int(0.5 + 10*(log($probmatrix[$c][$r] / $aminofreq[$c]))/log(10));
	$submatrix[$c][$r] = $submatrix[$r][$c] = int(0.5 + (10*(log($probmatrix[$c][$r] / $aminofreq[$c]))/log(10) + 10*(log($probmatrix[$r][$c] / $aminofreq[$r]))/log(10))/2);
    }
}

$outfilename = 'PAM' . $divergence . 'scores.txt';
open(OFH, '>', $outfilename) or die $!;

# Print out the substitution matrix!
print OFH ">PAM$divergence \n";
#print OFH "@keys\n";
print OFH "A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V\n";
for ($i = 0; $i < 20; $i++) {
#    print OFH "$keys[$i]: ";
    for ($j = 0; $j < 20; $j++) {
	if ($j < 19) {print OFH "$submatrix[$i][$j],";}
	else {print OFH "$submatrix[$i][$j]";}
    }
    print OFH "\n";
}
    


