#!/usr/local/bin/perl
BEGIN {
    print "Content-type: text/html\n\n";
    # Ensure that errors go to the web browser. From Steve's Primer3 program
    open(STDERR, ">&STDOUT");
    $| = 1;
    print '';
}

# Staffa path for this file is:
# /usr/local/etc/httpd/cgi-bin/pagelab/charles/

strict;
print <<EOF;
<html><TITLE>Fasta Data Modifier</TITLE>
<BODY BGCOLOR=lightblue  text="#000000" link="#118911">
EOF

require "ParseQuery.pl";
&PARSE;
&DEFINECONSTANTS;
&OUTPREP;
&LOADFASTA;
&PREPSUBMISSION;
	

# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub OUTADD {
	my $sts		= shift;	# [0] Name of STS
	my $left	= shift;	# [1] Left Primer
	my $right	= shift;	# [2] Right Primer
	my $base	= shift;	# [3] Base Sequence
	my $noBase	= shift;	# [4] Flag to indicate that consensus base sequence does not match.
	
	open (FILE, ">>$HOME/Data/SubmissionOut") or die "Failure to append to $HOME/Data/SubmissionOut: $!\n";
	
	print FILE "TYPE: STS\n";
	print FILE "STATUS: New\n";
	print FILE "CONT_NAME: Charles A. Tilford\n";
	print FILE "PROTOCOL: STS-RH (C. Tilford)\n";
	print FILE "BUFFER: Promega (CAT)\n";
	print FILE "CITATION:\n";
	print FILE "The Human Y Chromosome: A Fusion of Gene-laden Amplicons and X-homologous Domains\n";
	print FILE "SOURCE: YAC Subtraction\n";
	print FILE "STS#: $alias{$sts}\n";
	print "<FONT COLOR=red>No alias found for $sts!</FONT>\n" if ($alias{$sts} eq "");
	print FILE "SYNONYMS: $sts\n";
	print FILE "F_PRIMER: $left\n";
	print FILE "B_PRIMER: $right\n";
	print FILE "DNA_TYPE: Genomic\n";
	print FILE "PUBLIC: 08/01/2001\n";
	
	my $t = $temp{$sts};
	if ($t < 50) {
		print "<FONT COLOR=red>Error - $sts has a temperature of $t</FONT>\n";
		$t = 62;
	}
	print FILE "PCR_PROFILE:\n";
	print FILE "     94C        4:00 min\n";
	print FILE "         / 94C   :30 sec\n";
	print FILE "     35x | ${t}C  2:00 min\n";
	print FILE "         \\ 72C  2:00 min\n";
	print FILE "     72C        7:00 min\n";
	print FILE "SEQUENCE:\n";
	my $len = index($base, &REVERSE($right)) - index($base, $left);
	if ($noBase) {
		my $sp = substr($sts,0,7);
		print FILE &BLOCKOUT($Original{$sp}, 60);
		print FILE "COMMENT:\n";
		print FILE "The base sequence listed above was later found to contain\n";
		print FILE "errors. Although these primers still constitute a functional\n";
		print FILE "STS, we believe at this time that the following sequence more\n";
		print FILE "accurately reflects the base sequnce:\n";
		$len = index($Original{$sp}, &REVERSE($right)) - index($Original{$sp}, $left);
		print "<FONT COLOR=darkgray>$sts has been noted as having erroneous base sequence.</FONT>\n";
	}
	$len += length($right);
	print FILE &BLOCKOUT($base, 60);
	print FILE "SIZE: $len\n";
	print "<FONT COLOR=brick>$sts has a length of $len.</FONT>\n" if ($len < 75);
	print FILE "||\n\n";
	close FILE;
	$written{$sts}++;
	print "<FONT COLOR=red>WARNING!! $sts has $written{$sts}++ entries!\n</FONT>" if ($written{$sts}++ > 1);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub BLOCKOUT {
	my $stuff = shift;
	my $size	= shift;
	my $out = "";
	my $counter = 0;
	while ($counter < length($stuff)) {
		$out .= substr($stuff, $counter, $size) . "\n";
		$counter += $size;
	}
	return $out;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub OUTPREP {
	my $general = "";
	open (FILE, "$HOME/Data/SubmissionGeneral") or die "Failure to open $HOME/Data/SubmissionGeneral: $!\n";
	while (<FILE>) {
		$general .= $_;
	#	print "$_<BR>";
	}
	close FILE;

	open (FILE, ">$HOME/Data/SubmissionOut") or die "Failure to write to $HOME/Data/SubmissionOut: $!\n";
	print FILE $general;
	close FILE;
	
	print "<FONT COLOR=darkgreen>Finished reading base submission data.</FONT><BR>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Scours a file of FASTA sequences
sub LOADFASTA {
	@theNames = @theSeqs = ();
#	my $target = "OldFasta";		# File containing raw FastA entries
	my $target = "NextFasta";		# File containing sequencher-processed entries
	
	
	open (FILE, "$HOME/Data/$target") or die "Failure to open $HOME/Data/$target: $!\n";
	my $data = ""; my $name = "";
	while (<FILE>) {
		chomp;
		my @lin = split /[\n\r]/, $_;
		foreach my $l (@lin) {
			if ($l =~ />/) {
				$_ = $l;
				s/>//; s/<//;
				my $newName = $_;
				&STOREDATA($data, $name);
#				&TWIDLEDATA($data, $name);
				$data = "";
				$name = $newName;
			} else {
				$data .= $_;
			}
		}
	}
#	&TWIDLEDATA($data, $name);
	&STOREDATA($data, $name);
	close FILE;
#	&SUBFASTA;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PREPSUBMISSION {

	%left = %right = ();		#left and right primers
	open (FILE, "$HOME/Data/STSdata") or die "Failure to open $HOME/Data/STSdata: $!\n";
	while (<FILE>) {
		chomp;
		my @fields = split "\t", $_;
		my $sts = $fields[0];
		next if ($Ordered{$sts} eq "");	# Not a marker we're interested in
		print "<H1>$sts</H1>" if ($sts =~ /SP43-25/i);
		$left{$sts}		= $fields[3];
		$right{$sts}	= $fields[4];
	}
	close FILE;

	print "<PRE>";
	my @checking = sort {$a cmp $b} keys %left;
	foreach my $k (@checking) {
		my $sp = substr($k,0,7);
		if ($#{$location{$sp}} < 0) {
			print "<FONT COLOR=brown>No base sequences found for $sp:\n</FONT>";
			&FINDMATCHES(\@allIndices, $k, 1);					
		} else {
			&FINDMATCHES(\@{$location{$sp}}, $k);
		}
	}
	foreach my $k (keys %baseHits) {
		if ($#{$baseHits{$k}} > 0) {
			print "<FONT COLOR=purple>Multiple hits on the base sequence $k:</FONT>\n";
			foreach my $s (@{$baseHits{$k}}) {
				print "   $s\n";
			}
		}
	}
	print "Total fruity matches: $fruityMatch\n";
	print "Total lousy matches: $lousyMatch\n";
	
	print "</PRE>";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub FINDMATCHES {
	my $seqsRef	= shift;	# [0] - reference to list of indices for theSeqs
	my $sts		= shift;	# [1] - name of STS being checked
	my $verbose	= shift;	# [2] - indicates if total output should be reported
	foreach my $n (@{$seqsRef}) {
		my $lOrient = $rOrient = 1;
		my $lMatch = index($theSeqs[$n], $left{$sts});
		my $rMatch = index($theSeqs[$n], $right{$sts});
		if ($lMatch < 0) {
			$lMatch = index($theSeqs[$n], &REVERSE($left{$sts}));
			$lOrient = -1;
			$lOrient = 0 if ($lMatch < 0);
		}
		if ($rMatch < 0) {
			$rMatch = index($theSeqs[$n], &REVERSE($right{$sts}));
			$rOrient = -1;
			$rOrient = 0 if ($rMatch < 0);
		}
		my $prod = $lOrient * $rOrient;
		my $for = $left{$sts}; my $rev = $right{$sts};
		if ($lOrient == -1 || $rOrient == 1) {
			$rev = $left{$sts}; $for = $right{$sts};
		}
		if ($prod == -1) {		# Both primers match, opposite orientations
			if ($lMatch * $lOrient < -1 * $rMatch * $rOrient) {
				print "<FONT COLOR=green>$sts: Solid double hit on $theNames[$n].\n</FONT>" if ($verbose);
				push @{$baseHits{$theNames[$n]}}, "$sts Both";
				&OUTADD($sts, $for, $rev, $theSeqs[$n]);
			} else {
				print "$sts: Both primers match $theNames[$n], <FONT COLOR=red>but the orientation is incorrect.</FONT>\n";
				push @{$baseHits{$theNames[$n]}}, "$sts Both Backwards";
			}
			last;
		} else {
			if ($lOrient != 0) {
				my $checking = $right{$sts}; my $Or = 1;
				if ($lOrient == 1) {		# Make sure to look in the opposite direction
					$checking = &REVERSE($right{$sts});
					$Or = -1;
				}
				my @Positions = ();
				my $misMatch = &NEARMATCH(\$checking, \$theSeqs[$n], \@Positions);
				if ($#Positions == 0) {
		#			print "$sts: <FONT COLOR=blue>Good left match to $theNames[$n], </FONT> Poor match with $misMatch mismatch at $Positions[0] bp:\n";
		#			&SHOWMIS($left{$sts}, $lMatch, $lOrient, $right{$sts}, $Positions[0], $Or, $theSeqs[$n]);
					push @{$baseHits{$theNames[$n]}}, "$sts Left Right=$misMatch";
					&OUTADD($sts, $left{$sts}, $right{$sts}, $theSeqs[$n], 1);
					$fruityMatch++;
				} else {
					print "$sts: <FONT COLOR=blue>Good left match to $theNames[$n], </FONT>";
					printf ("<FONT COLOR=red>%d mismatches on right primer, in %d positions. </FONT>\n", $misMatch, $#Positions + 1);
					$lousyMatch++;
				}
			} elsif ($rOrient != 0) {
				my $checking = $left{$sts}; my $Or = 1;
				if ($rOrient == 1) {
					$checking = &REVERSE($left{$sts});
					$Or = -1;
				}
				my @Positions = ();
				my $misMatch = &NEARMATCH(\$checking, \$theSeqs[$n], \@Positions);
				if ($#Positions == 0) {
		#			print "$sts: <FONT COLOR=blue>Single right match to $theNames[$n].</FONT> Poor match with $misMatch mismatch at $Positions[0] bp:\n";
		#			&SHOWMIS($left{$sts}, $Positions[0], $Or, $right{$sts}, $rMatch, $rOrient, $theSeqs[$n]);
					push @{$baseHits{$theNames[$n]}}, "$sts Right Left=$misMatch";
					&OUTADD($sts, $left{$sts}, $right{$sts}, $theSeqs[$n], 1);
					$fruityMatch++;
				} else {
					print "$sts: <FONT COLOR=blue>Single right match to $theNames[$n]. </FONT>";
					printf ("<FONT COLOR=red>%d mismatches on left primer, in %d positions. </FONT>\n", $misMatch, $#Positions + 1);
					$lousyMatch++;
				}
			} else {
				# print "<FONT COLOR=red>Screw Up! No matches for $theNames[$n].\n</FONT>";
			}
		}
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub SHOWMIS {
	my $lPri	= shift;
	my $lPos	= shift;
	my $lOr		= shift;
	my $rPri	= shift;
	my $rPos	= shift;
	my $rOr		= shift;
	my $base	= shift;
	$lPri		= &REVERSE($lPri) if ($lOr == -1);
	$rPri		= &REVERSE($rPri) if ($rOr == -1);
	
	my $leftmost = $lPos;
	$leftmost = $rPos if ($rPos < $lPos);
	if ($leftmost > 5) {
		my $delta = $leftmost - 5;
		$lPos += 3 - $delta;
		$rPos += 3 - $delta;
		$base = "..." . substr($base, $delta);
	}
	my $endPoint = $rPos + length($rPri) + 5;
	$endPoint = $lPos + length($lPri) + 5 if ($rPos < $lPos);
	$base = substr($base, 0, $endPoint) . "..." if (length($base) > $endPoint);
	
	my $left	= "<FONT COLOR=green>";
	$left		.= " " x $lPos;
	my $l = length($lPri) - 1;
	for my $a (0..$l) {
		my $col = "";
		if (substr($base,$lPos+$a,1) ne substr($lPri,$a,1)) {
			$col = "<FONT COLOR=red>";
			$col = "<FONT COLOR=purple>" if (substr($lPri,$a,1) =~ /[$Ambig{substr($base,$lPos+$a,1)}]/);
		}
		$left .=  $col . substr($lPri,$a,1);
		$left .=  "</FONT>" if ($col);
	}
	$left .= "</FONT>";
	my $right	= "<FONT COLOR=green>";
	$right		.= " " x $rPos;
	my $r = length($rPri) - 1;
	for my $a (0..$r) {
		my $col = "";
		if (substr($base,$rPos+$a,1) ne substr($rPri,$a,1)) {
			$col = "<FONT COLOR=red>";
			$col = "<FONT COLOR=purple>" if (substr($rPri,$a,1) =~ /[$Ambig{substr($base,$rPos+$a,1)}]/);
		}
		$right .=  $col . substr($rPri,$a,1);
		$right .=  "</FONT>" if ($col);
	}
	$right .= "</FONT>";
	print "$left\n";
	print "<FONT COLOR=BLACK>$base</FONT>\n";
	print "$right\n\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub NEARMATCH {
	my $priRef	= shift;	# [0] reference to string holding primer sequence
	my $baseRef	= shift;	# [1] reference to base sequence string
	my $locRef	= shift;	# [2] reference to array holding results
	my $priLength = length($$priRef) - 1;
	my $last = length($$baseRef) - $priLength - 2;
	my @hits = ();
	for my $b (0..$last) {
		my $score = 0;
		for my $p (0..$priLength) {
			$score++ if (substr($$baseRef, $b+$p, 1) ne substr($$priRef, $p, 1));
		}
		push @hits, sprintf("%4.4s\t%d", $score, $b);
	}
	@hits = sort {$a cmp $b} @hits;
	my ($bestScore, $bestPos) = split "\t", $hits[0];
	@{$locRef} = ($bestPos);
	my $counter = 0; my ($nextS, $nextP);
	do {
		$counter++;
		($nextS, $nextP) = split "\t", $hits[$counter];
		push @{$locRef}, $nextP if ($nextS == $bestScore);
	} while ($nextS == $bestScore);
	return $bestScore + 0;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub REVERSE {
	my $d = $_[0];
	my $rev = "";
	while (length($d) > 0) {
		$n = chop $d;
		if ($rc{$n}) {
			$rev .= $rc{$n};
		} else {
			$rev .= "#";
		}
	}
	return $rev;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub STOREDATA {
	my $seq			= shift;
	my $theName		= shift;
	push @theNames, $theName;
	push @theSeqs, $seq;
	push @allIndices, $#theNames;
	my @SPs = split '/', $theName;
	for my $s (0..$#SPs) {
		my $sp = $SPs[$s];
		$sp = "SP" . $sp if ($s > 0);
		$sp = substr($sp,0,7);
		push @{$location{$sp}}, $#theNames;
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Generate a bunch of potentially smaller FASTA entries, split at a defined locus and semi-filtered
sub SUBFASTA {
	my $num = 150;
	my $count = 0;
	while ($count * $num < $#theNames) {
		my $tag = $count;
		$tag = "0" . $tag if ($tag < 10);
		my $filename = "NewFasta" . $tag;
		open (FILE, ">$HOME/Data/FastaOut/$filename") or die "Failure to write to $HOME/Data/FastaOut/$filename: $!\n";
		$low = $count * $num;
		$high = $low + $num - 1;
		$high = $#theNames if ($high > $#theNames);
		for my $i ($low..$high) {
			print FILE ">$theNames[$i]\n";
			print FILE "$theSeqs[$i]\n";
		}
		close FILE;
		system "chmod 666 $HOME/Data/FastaOut/$filename";
		$count++
	}
	print "Total entries = $#theNames<BR>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# This routine was only used early on to automatically chop sequences up at likely chimeric sites.
sub TWIDLEDATA {
	return if ($_[0] eq "");
	$_ = $_[0];
	foreach my $in (@swapping) {
		$out = "$in\n$in";
		s/$in/$out/g;
	}
#	print "$_[1] $_<BR>";
	my @lines = split "\n", $_;
	my $tag = 65;
	my @keep = ();
	foreach my $l (@lines) {
		my $len = length($l);
		next if ($len < 30);
		my @tmpA = ($l =~ /(N)/g);
		$totN = $#tmpA + 1;
		next if ($totN/$len > 0.25);
		push @keep, $l;
	}
	foreach my $l (@keep) {
		my $theName = $_[1];
		$theName .= "-" . chr($tag) if ($#keep > 0);
		$theName .= "*" if ($Ordered{$_[1]} ne "");
		push @theNames, $theName;
		$tag++;
		push @theSeqs, $l;
		
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub DEFINECONSTANTS {
	$HOME = "/usr/local/etc/httpd/cgi-bin/pagelab/charles";
	$Here = "."."$ENV{'REMOTE_ADDR'}";				# Returns the IP address of the *USER'S* Machine (to allow multiple users)
	$DataLoc =	"/tmp/ChasFasta" . $Here;			# Unique directory for the user in the tmp folder
	$WebLoc =	"/tmp/ChasFasta" . $Here;			# The base URL for the host computer
	@swapping = ("GATC", "NATC", "GNTC", "GANC", "GATN");
	
	$Allowed	= "AGCTUYRWSKMBDHVXN";
	$Reverse	= "TCGAARYWSMKVHDBXN";
	%Ambig 		= ( A => 'A',	M => 'AC',	S => 'CG',	V => 'ACG',	X => 'ACGT',
					C => 'C',	R => 'AG',	Y => 'CT',	H => 'ACT',	N => 'ACGT',
					T => 'T',	W => 'AT',	K => 'GT',	D => 'AGT',
					G => 'G',							B => 'CGT', 				);
	%rc = ();
	for $i (0..length($Allowed)-1) {
		$rc{ substr( $Allowed , $i ,1) } = substr( $Reverse , $i ,1 );
	}
	
	%Ordered = ();		# STSs that have been ordered (ie purchased)
	open (FILE, "$HOME/Data/NotOnY") or die "Failure to open $HOME/Data/NotOnY: $!\n";
	while (<FILE>) {
		chomp;
		my ($name, $location) = split "\t", $_;
		$NotOnY{$name} = "Y";
	}
	close FILE;
	
	open (FILE, "$HOME/Data/SPaliases") or die "Failure to open $HOME/Data/SPaliases: $!\n";
	while (<FILE>) {
		chomp;
		my ($sts, $sY) = split "\t", $_;
		$alias{$sts} = $sY;
	}
	close FILE;
	
	open (FILE, "$HOME/Data/IdenticalSP") or die "Failure to open $HOME/Data/IdenticalSP: $!\n";
	while (<FILE>) {
		chomp;
		my ($discarding, $keeping) = split "\t", $_;
		$Identical{$discarding} = $keeping;
	}
	close FILE;
	
	open (FILE, "$HOME/Data/OriginalSP") or die "Failure to open $HOME/Data/OriginalSP: $!\n";
	while (<FILE>) {
		chomp; s/>//;
		my $sts = $_;
		my $seq = <FILE>; chomp $seq;
		$Original{$sts} = $seq;
	}
	close FILE;
	
	@qq = ("Excellent", "Reasonable", "Poor", "Useless", "");
	for my $i (0..$#qq) {
		$qualVal{$qq[$i]} = $i;
	}
	%temp = ();
#	print "<PRE>";
	open (FILE, "$HOME/Data/AllQuality") or die "Failure to open $HOME/Data/AllQuality: $!\n";
	while (<FILE>) {
		chomp;
		my ($sts, $topScore, $topTemp, $botScore, $botTemp) = split "\t", $_;
		$temp{$sts} = $topTemp;
		$temp{$sts} = $botTemp if ($qualVal{$botScore} < $qualVal{$topScore});
#		print "$sts\t$temp{$sts}\n" if ($sts =~ "sY");
	}
	close FILE;
	
	open (FILE, "$HOME/Data/OrderedSP") or die "Failure to open $HOME/Data/OrderedSP: $!\n";
	while (<FILE>) {
		chomp;
		my ($name, $alias) = split "\t", $_;
	#	print "$name<BR>" if ($name =~ /[ABC]/);
		next if (exists $NotOnY{$name});
		next if (exists $Identical{$name});
	#	$name =~ s/\D$//;
		next if (exists $NotOnY{$name});
		next if (exists $Identical{$name});
		$Ordered{$name} = "Y";
	}
	close FILE;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
