#jacob hanger

&partA;
&partB;



sub partA  #this part opens the newGenome and  file and determines the relative frequencies of each amino acid in each of th
{
	#	a	t	c	g
#START ALL THE VARAIBLES THAT I WILL USE IN THE REST OF THE PROGRAM
	($CGrichsequence, $backgroundsequence) = &getFiles;
	$CGrichsequence = lc($CGrichsequence);
	$backgroundsequence = lc($backgroundsequence);
	$CGrichrichtally = 	$backgroundtally = 0;
	$tnstance = $CGrichsequence;
	$tnstance =~ s[(.)(?=.*?\1)][]g; #reduces the string to only the composing characters "aaac"=> "ac"
 	@character = split(//, $tnstance);
	@CGrichseq = split(//, $CGrichsequence);
	@backgroundseq = split(//, $backgroundsequence);
	
	
	
	for $t (0 ... (scalar(@character)-1))	 #initialize EVERYTHING.
	 {
			$backgroundprob[$t]  = $CGrichprob[$t] = 0;	
	}
	
	for($t = 0; $t<scalar(@CGrichseq)-1; $t++) #tallies up the entire sequences
	{
		$CGrichx = &indx($CGrichseq[$t]);
		$CGrichrichtally++;
		$CGrichprob[$CGrichx]+=1;
	}
	
	for($t = 0; $t<scalar(@backgroundseq)-1; $t++) #tallies up the entire sequences
	{
		$backgroundx = &indx($backgroundseq[$t]);
		$backgroundtally++;
		$backgroundprob[$backgroundx]+=1;
	}

	
	for $x (0 ... (scalar(@character)-1))	#goes thru and makes it a percentage. if the number is 0, it makes it a very small number
	{
			$CGrichprob[$x] = $CGrichprob[$x]/$CGrichrichtally;
			$backgroundprob[$x] = $backgroundprob[$x]/$backgroundtally;
			if($CGrichprob[$x] eq 0)
			{
				$CGrichprob[$x] = exp(-145);
			}
			if($backgroundprob[$x] eq 0)
			{
				$backgroundprob[$x] = exp(-145);
			print "$backgroundprob[$x]      ";
			}
	}
	print "\n";
}


sub partB() #this is the same method as I use in problem 2. 
{
@data = split(//, lc(&initialize));
%probdeEMIS = ("C" => { "t" => $CGrichprob[1], "a" => $CGrichprob[0], "c" => $CGrichprob[2], "g" => $CGrichprob[3] },  "B" =>{ "t" => $backgroundprob[1],"c" => $backgroundprob[2],"a" => $backgroundprob[0],  "g" => $backgroundprob[3] });
@estados=("C", "B");
%startProbOfestado=("C" => .5 , "B" => .5);
%probdeCam = ("C" => {"C" => .99, "B" => .01}, "B"=> {"C" => .001, "B" => .999});


foreach(@estados)
{	
	@database[$tndx]=[$startProbOfestado{$_}, $_, $startProbOfestado{$_}];
	$tndx++;
}

$tndx=0;
foreach(@data)
{
   	for ($siguiente_subestado=0; $siguiente_subestado < @estados; $siguiente_subestado++)
   	{
       	$mayor=-exp(99999);
       	for ($subestado=0; $subestado<@estados; $subestado++)
		{	
			$prob=$database[$subestado][0];
			$v_path=$database[$subestado][1];
			$vertaaaahaa=$database[$subestado][2];
			$p1= log($probdeEMIS{$estados[$subestado]}{$_});
			$p2 = log($probdeCam{$estados[$subestado]}{$estados[$siguiente_subestado]});
			$p = $p1+$p2;
           	$prob=$prob + $p;
           	$vertaaaahaa=$vertaaaahaa + $p;
           	$total=$total + $prob;
           	if ($vertaaaahaa > $mayor)
           	{   
				$listamax = $v_path."".$estados[$siguiente_subestado];
				$mayor=$vertaaaahaa;
			}
		}
		@snapshot[$siguiente_subestado]= [$total, $listamax, $mayor];
	}
	@database=@snapshot;
}

$listamax=$database[0][1];
$mayor=$database[0][2];
for ($t=0; $t<@estados; $t++)
{
   	$prob=$database[$t][0];
   	$v_path=$database[$t][1];
   	$vertaaaahaa=$database[$t][2];
   	$total=$total + $prob;
   	if ($vertaaaahaa > $mayor)
      	{
          	$listamax=$v_path;
          	$mayor=$vertaaaahaa
      	}
}

@data = split(//, $listamax);
@indices;
for($t = 0; $t<scalar(@data)-1; $t++)
{
	$actual = $data[$t];
	$siguiente = $data[$t+1];
	if($actual eq 'B' and $siguiente eq 'C')
	{
		push(@indices, $t);
	}
}

print "The CG-rich regions begin at ";
print @indices;
print ".\n\n"
}

sub initialize()
{
	open FILE, "newGenome.fasta" or die("Could not open file!");
	@data=<FILE>;
	$string ="";
	$pila = 0;
	$probdeCamCC = .99;
	$probdeCamCB = .01;
	$probdeCamBB = .999;
	$probdeCamBB = .001;
	$startProbC = $startProbB = .5;
	foreach(@data)
	{
		chomp;
		if($pila eq 0)
		{
			$pila++;
		}
		else
		{
		$string .= $_;
		}
	}
	return $string;	
	
}

sub getIndex
{
	@y = @_[0];
	@indices;
	for($t = 0; $t<scalar(@y)-1; $t++)
	{
		$actual = $y[$t];
		$siguiente = $y[$t+1];
		if($actual eq 'B' and $siguiente eq 'C')
		{
			push(@indices, $t);
		}
	}
	
}

sub getFiles
{
	open CGFILE, "newGenome.fasta" or die("Could not open file!");
	my	$CG ="";
	my	$BG ="";
	my	@CGrichData=<CGFILE>;
	my	$pila = 0;
	foreach(@CGrichData)
	{
		chomp;
		$st = $_;
		if($pila >0)
		{
		$CG .= $_;
		}
	$pila++;
	}
	open backgroundFILE, "newBackgroundSequences.fasta" or die("Could not open file!");
	my	@backgroundData=<backgroundFILE>;
	foreach(@backgroundData)
	{
		chomp;
		$BG .= $_;	
	}
	
	return ($CG, $BG);
}

sub indx()
{
	my $backgroundx = 5;
	my $backgroundactual = @_[0];
	if($backgroundactual eq 'a')
	{
		$backgroundx=0;
	}
	if($backgroundactual eq 't')
	{
		$backgroundx=1;
	}
	if($backgroundactual eq 'c')
	{
		$backgroundx=2;
	}
	if($backgroundactual eq 'g')
	{
		$backgroundx=3;
	}
	return $backgroundx;
}