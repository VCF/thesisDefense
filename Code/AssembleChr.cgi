#!/usr/local/bin/perl
BEGIN {
    print "Content-type: text/html\n\n";
    # Ensure that errors go to the web browser. From Steve's Primer3 program
    open(STDERR, ">&STDOUT");
    push @INC, "/usr/local/etc/httpd/cgi-bin/pagelab/charles/FindOligos";
    $| = 1;
    print '';
}


# Staffa path for this file is:
# /usr/local/etc/httpd/cgi-bin/pagelab/charles/
################################################
# Required Directory Structure and Support Files
# - ParseQuery.pl		Routine for parsing HTML arguments
# - CGI Package			Lincoln Stein's CGI library
# - GD Package			Lincoln's GIF graphics manipulation library
# Data/					Folder at same level as Fasta.cgi
# Data/RepeatLibrary	Fasta Library of genomic repeats

use CGI;
use CGI::Carp 'fatalsToBrowser';
use GD;
use FindOligos;
require "ParseQuery.pl";
require "ParseFasta.pl";
srand ( time() ^ ($$ + ($$<<15)) );								# Make a "random" seed for random number generation


&DEFINECONSTANTS;
&PARSE;

print "<!\n";													# Comment out any annoying messages about rm
unless ($AlignThem ne "" || $File ne "") {						# Don't delete the temporary files if the program was run as a pop-up alignment report
	system "rm $DataLoc/*.tmp";									# Remove old gif files
	$theCommand = "rm $DataLoc/*.aln";
	print "\nExecuting $theCommand: "; system "$theCommand";	# Occasionally there will be too many files to delete - the system will choke and just report the arg list being too long
}
print ">\n";

print "<html><TITLE>Genomic Clone Assembler</TITLE><BODY BGCOLOR=White>";
#&PRINTOPTIONS;

# my @variants = ();
# &MAKEAMBIG("ATCGATCGTTACG", 3, \@variants);

&SPLICEBACS if ($splice);
&PLACESTS if ($place);
&DRAWOPTIONS;

# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub DRAWOPTIONS {
	print "<FONT COLOR=orange SIZE=+2>Genomic Assembler Options</FONT><BR>\n";
	print qq(<P><A HREF="AssembleChr.cgi?splice=Y&rand=$$">Splice Bacs</A> - Find overlaps in BAC collection and splice them into larger contigs.\n);
	print qq(<P><A HREF="AssembleChr.cgi?place=Y&rand=$$">ePCR</A> - Electronically identify STS sites.\n);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PLACESTS {
	my $dir = "$HOME/Data/ChrYChunks";		# The directory where the clones lurk
	my @temp = split "\n", `ls $dir`;
	my @fileList = ();
	my @toSkip = ( );
	foreach $s (@toSkip) { $skipping{$s} = "Y";}
	foreach $f (@temp) {
		next if (exists $skipping{$f});				# Redundant
		push @fileList, $f;
	}
	$| = 1;
	print "<CENTER><TABLE BGCOLOR=yellow BORDER=1><TR WIDTH=300><TD><CENTER><H3>STS Identification by ePCR</H3>\n";
	print "<FONT COLOR=red><H2>Long Calculation in Progress<BR>Please do not Disturb</H2></FONT></CENTER></TD></TR></TABLE></CENTER>\n";


#	open (FILE, "$HOME/Data/PrimerData") or die "PLACESTS: Can't read PrimerData: $!\n";
#	undef (%left); undef (%right); undef (%size);
#	while (<FILE>) {
#		chomp;
#		my ($sts, $left, $right, $size) = split "\t", $_;		
#		$left{$sts} = uc($left);
#		$right{$sts} = uc($right);
#		$size{$sts} = $size;
#	}
#	close FILE;

	undef (%left); undef (%right); undef (%size);
	open (FILE, "$HOME/Data/_STSData") or die "Failure to open _STSData: $!\n";
	while (<FILE>) {
		chomp;
		my ($sts, $d1, $d2, $left, $right, $d3, $d4, $size) = split "\t", $_;
		$sts =~ s/SQ/SP/;
		$name = $sts; $name = $alias{$sts} if ($alias{$sts} ne "");
		$left{$name} = uc($left);
		$right{$name} = uc($right);
	}
	close FILE;
	
	undef (%accession);
	open (FILE, "$HOME/Data/HelenSTSdata") or die "Failure to open HelenSTSdata: $!\n";
	my $sts = "";
	while (<FILE>) {
		chomp;
		($key,$val) = split " ", $_;
		next if ($val eq "");
		$sts = $val if ($key =~ /STS_name/);
		next if ($sts eq "");
		$left{$sts} =  uc($val)	if ($key =~ /Sequence_left/);
		$right{$sts} = uc($val)	if ($key =~ /Sequence_r/);
		$size{$sts} = $val		if ($key =~ /size/i);
	}
	close FILE;

	$lowProdSize = 50;
	$highProdSize = 1000;
	
	foreach $f (@fileList) {
		print "<FONT COLOR=green SIZE=+2><B>$f</B></FONT><BR>\n";
		my @Sequences = my @Names = ();
		my $numLines = &LOADFASTA("$dir/$f", \@Sequences, \@Names);
		my $name = $Names[0];
		my $seq = $Sequences[0];
		my $length = length ($seq);
		my %Matches = ();
		my $markersChecked = 0;

		FindOligos::setTemplate($seq);
		my $t = time;
		$ambigs = 1;
		foreach $sts (keys %left) {
			next if ($left{$sts} eq "" or $right{$sts} eq "");
			@{$Matches{$sts}} = FindOligos::findSts($left{$sts}, $right{$sts}, $ambigs, $lowProdSize, $highProdSize);
			$sites{$sts} = shift @{$Matches{$sts}};
	#		$sites{$sts} = &FINDSTS(\$seq, $left{$sts}, $right{$sts}, $size{$sts}, 0, \@{$Matches{$sts}}, $lowProdSize, $highProdSize);
			$markersChecked++;
			print " $sts=$sites{$sts}";				# Keep the browser interested...
	#		last if ($markersChecked > 4);
		}
		print "<BR>\n";
		my @s = ();
		foreach $sts (keys %Matches) {
			push @s, sprintf("%9s\t%s", $Matches{$sts}[0], $sts) unless ($Matches{$sts}[0] eq "");
		}
		@s = sort {$a cmp $b} @s;
		
		my $fName = "${f}Markers-test";
		open (FILE, ">$HOME/Data/ChrYMarkers/$fName") or die "LOADFILES: Can't write to $fName: $!\n";
		print FILE "$f\t$length\n";
		foreach $i (@s) {
			my ($dumb, $sts) = split "\t", $i;
			print "<B>$sts</B>\t<i>$sites{$sts}</i>";
			print FILE "$sts\t$sites{$sts}";
			foreach $pos (@{$Matches{$sts}}) {
				print "\t$pos";
				print FILE "\t$pos";
			}
			print "<BR>\n";
			print FILE "\n";
		}
		$| = 1;
		close FILE;
		$t = time - $t;
		my $rate = 0;
		if ($t > 0) {
			$rate = $markersChecked * $length/(1000 * $t);
		}
		printf ("<FONT COLOR=brown>%.0f seconds for %.0f markers and %.0f kb = %.2f kb-markers per second</FONT><BR>\n", $t, $markersChecked, $length/1000, $rate);
#		die;
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub FINDSTS {
	my $targetRef	= shift;		# [0] Reference to target sequence
	my $leftPri		= shift;		# [1] Left (forward) primer sequence
	my $rightPri	= shift;		# [2] Right (reverse) primer sequence
	my $prodSize	= shift;		# [3] Expected product size of primer
	my $ambiguous	= shift;		# [4] Number of mismatches to allow
	my $outputRef	= shift;		# [5] Array to hold matches for this STS
	my $lowSize		= shift;		# [6] lowest allowed product size
	my $hiSize		= shift;		# [7] Largest allowed product size
	$rightPri = &REVERSE($rightPri);
	
	$targLen = length($$targetRef);
	
	# Find oligo binding sites in both frames
	my @oligos	= ();		# [0] = forward frame, [1] = reverse frame
	my @Matches	= ();		# Will store base pair positions of matched primers
	push @oligos, [( $leftPri,	&REVERSE($rightPri) )];	# Forward frame primers
	push @oligos, [( $rightPri,	&REVERSE($leftPri) )];	# Reverse frame primers

    $| = 1;

	for $direction (0..1) {								# Check both frames
		foreach $primer (@{$oligos[$direction]}) {		# Check both primers in the frame
			my $priLen = length($primer);
			my $t = time;
			if ($ambiguous == 0) {
				my $zeroPoint = $[; $zeroPoint++;	# ] just a bracket to help BBEdit balance parentheses
				my $looking = $zeroPoint;
				while (($looking = index($$targetRef, $primer, $looking+1)) >= $zeroPoint) {
					push @{$Matches[$direction]}, sprintf("%7s", $looking + $priLen*$direction);
				}
			} else {
				for $targPos (0..($targLen-$priLen)) {		# Scan along entire target length
					my $misMatch = $priPos = 0;
					my $busyWork = substr($$targetRef, $priPos+$targPos, $priLen);
					next;
					while ($misMatch <= $ambiguous && $priPos < $priLen) {
						$misMatch++ if (substr($primer, $priPos, 1) ne substr($$targetRef, $priPos+$targPos, 1));
						$priPos++;
					}
					if ($misMatch <= $ambiguous) {
						push @{$Matches[$direction]}, sprintf("%7s", $targPos + $priLen*$direction);
						# The position is increased by primer length in the reverse frame - standardize to 5' nucleotide
						# Will be off by 1 for reverse frame
					}
				}
			}
		}
		@{$Matches[$direction]} = sort {$a cmp $b} @{$Matches[$direction]};
	}
	# Find oligo pairs that would be expected to produce PCR products - we also allow an olgio to amplify off itself
	push @{$Matches[1]}, 9999999999999;					# A big number to serve as stopper
	my $revIndex = my $size = 0;
	for $i (0..$#{$Matches[0]}) {						# Cycle through ordered forward frame
		while ($Matches[1][$revIndex] < $Matches[0][$i]) {						# Time saver (?) - only check primers to "right" of forward primer
			$revIndex++;
		}
		for $j ($revIndex..$#{$Matches[1]}) {
			last if (($size = $Matches[1][$j] - $Matches[0][$i]) > $hiSize);	# Stop checking if the reverse primers are too far away
			push @{$outputRef}, ($Matches[0][$i] + int($size/2)) unless ($size < $lowSize);
			# The position of the STS is taken to be the center of the amplified region
		}
	}
	return ($#{$Matches[0]} + $#{$Matches[1]} + 1);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Takes a sequence and generates all possible sequences with $variants ambiguities
sub MAKEAMBIG {
	my $baseSeq		= shift;
	my $variants	= shift;
	my $outputRef	= shift;
	
	
	my $len = length($baseSeq);
	my $lenIn = $len - 1;
	my @count = ();
	
	my $vars = $variants - 1;
	for $i (0..$vars) {
		$count[$i] = $i;						# These are the positions that are currently being "ambigified" - made "N" (or "." in match terminology)
	}
	$count[$vars+1] = $len*2;					# Dummy counter, with large value as "stopper"
	while ($count[$vars] < $len) {				# So long as the last spot has room to move right...
		while ($count[$vars] < $len) {		
			push @{$outputRef}, &VARIANT($baseSeq, \@count);
			$count[$vars]++;					# Slide the right-most counter along to the right
		}
		my $advance = 0;
		until ($count[$advance] + 1 < $count[$advance+1]) {	# Find the first counter that is not immediately adjacent to its neighbor
			$advance++;
		}
		$count[$advance]++;						# Move the first available counter forward
		$count[$vars] = $count[$vars-1] + 1;	# Reset right counter to just past its leftward neighbor
	}
	@{$outputRef} = ($baseSeq) if ($#{$outputRef} == -1);
	print "<PRE>";
	foreach $v (@{$outputRef}) {
		print "$v\n";
	}
	die;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub VARIANT {
	my $baseSeq		= shift;	# The basic sequence
	my $countRef	= shift;	# Array of the positions in above string to be replaced with "."
	for $i (0..($#{$countRef}-1)) {
		substr($baseSeq, ${$countRef}[$i], 1) = ".";
	}
	return $baseSeq;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub SPLICEBACS {
	print "<CENTER><FONT COLOR=green SIZE=+2>Splicing together BACs</FONT></CENTER></BR>\n";
	my $dir = "$HOME/Data/ChrYSequence";		# The directory where the clones lurk
	my @temp = split "\n", `ls $dir`;
	my @fileList = ();
	my @toSkip = ("0684N02", );
	foreach $s (@toSkip) { $skipping{$s} = "Y";}
	foreach $f (@temp) {
		next if (exists $skipping{$f});				# Redundant
		push @fileList, $f;
	}
	my $biteSize = 50;
	foreach $f (@fileList) {
		my @Sequences = my @Names = ();
		my $numLines = &LOADFASTA("$dir/$f", \@Sequences, \@Names);
		my $name = $Names[0];
		my $seq = $Sequences[0];
		my $length = length ($seq);
		$leftbit{$f}	= substr($seq,0,$biteSize);
		$rightbit{$f}	= substr($seq,$length-$biteSize);
	}
	print "<FONT COLOR=blue>Tailed Hashes constructed</FONT><BR>\n";
	
	open (FILE, "$HOME/Data/BACorder") or die "LOADFILES: Can't read BACorder: $!\n";
	undef (%Leftmost);
	while (<FILE>) {
		chomp;
		my @lines = split "\t", $_;
		$pos = shift @lines;
		foreach $bac (@lines) {
			$Leftmost{$bac} = $pos unless (exists $Leftmost{$bac});
		}
		if ($#lines >= 0) {
			push @totalLines, [@lines];
		}
	}
	close FILE;
	my $neighbors = 2;			# Number of neighboring markers that we are willing to look towards to hunt for overlaps
	for $i (0..$#totalLines) {
		my $l = $i - $neighbors;	$l = 0				if ($l < 0);
		my $r = $i + $neighbors;	$r = $#totalLines	if ($r > $#totalLines);
		my @localGroup = (); my %inthere = ();
		for $j ($l..$r) {
			foreach $bac (@{$totalLines[$j]}) {
				next unless (exists $leftbit{$bac});
				push @localGroup, $bac unless (exists $inthere{$bac});
				$inthere{$bac} = "Y";
			}
		}
		foreach $hip (@{$totalLines[$i]}) {
			next unless (exists $leftbit{$hip});
#			print "$hip, ";
			foreach $hop (@localGroup) {
				next if ($hip eq $hop);
				$neighbor{$hip}{$hop} = "Y";
			}
		}
#		print "<BR><BR>";
	}
	undef(@totalLines);
	
	print "<FONT COLOR=blue>Neighboring BACs identified</FONT><BR>\n";
	
	print "<PRE>";
	my $zeroPoint = $[;
	$zeroPoint++;
	print "Base index is set to $zeroPoint\n";
	foreach $f (@fileList) {
		print "$f\n";
		my @Sequences = my @Names = ();
		my $numLines = &LOADFASTA("$dir/$f", \@Sequences, \@Names);
		my $name = $Names[0];
		my $seq = $Sequences[0];
		my $length = length ($seq);
		my @locals = sort {$a cmp $b} keys %{$neighbor{$f}};
		foreach $k (@locals) {
			next if ($k eq $f);
			my $match = rindex($seq, $leftbit{$k});
			if ($match >= $zeroPoint) {
				print "   $k matches on RIGHT side at base $match\n";
				push @{$toRight{$f}}, $k;
			}
			$match = index($seq, $rightbit{$k});
			if ($match >= $zeroPoint) {
				print "   $k matches on LEFT  side at base $match\n";
				push @{$toLeft{$f}}, $k;
			}
		}
	}
	@{$toRight{"0055O11"}} = ("0538M13");	# Specify one of two options for overlapping
	@{$toRight{"0539D10"}} = ("0363G06");	# Palindromic area, false overlap with Bac to left
	@{$toRight{"0209I11"}} = ("0220O02");	# Overlap, but were too many markers away for the bacs to be considered by the algorithm
	@{$toLeft{"0220O02"}} = ("0209I11");	# Set Right match to above
	@{$toRight{"0477B05"}} = ("0427G18");	# Overlap, but were too many markers away for the bacs to be considered by the algorithm
	@{$toLeft{"0427G18"}} = ("0477B05");	# Set Right match to above
	@{$toRight{"0263A15"}} = ("0070G12");	# Overlap, but were too many markers away for the bacs to be considered by the algorithm
	@{$toLeft{"0070G12"}} = ("0263A15");	# Set Right match to above
	@{$toLeft{"0489O13"}} = ();				# This area is actually a gap, not a true overlap
	%ButtJoin = ();							# These are BACs that overlap only by the restriction site, but have been verified as overlapping from other sequence
	$ButtJoin{"0268K13"} = "0157F24";	
	$ButtJoin{"0251M08"} = "0100J21";	
	$ButtJoin{"0235I01"} = "0292E08";	
	
	foreach $left (keys %ButtJoin) {		# Set overlaps for butt joined BACs
		@{$toRight{$left}} = ($ButtJoin{$left});
		@{$toLeft{$ButtJoin{$left}}} = ($left);
	}
	@GrandContigs = ();
	foreach $f (@fileList) {		# Find starts of contigs then nucleate out from them
		next unless ($#{$toLeft{$f}} == -1 && $#{$toRight{$f}} <= 0);
		my @contig = ($f);
		my $focus = $f;
		while ($#{$toRight{$focus}} == 0) {
			$focus = $toRight{$focus}[0];
			push @contig, $focus;
		}
		push @GrandContigs, [@contig];
	}
	
	printf ("<FONT COLOR=green SIZE=+2>%d Total Contigs:</FONT>\n", $#GrandContigs + 1);
	undef(@Ordered);
	for $i (0..$#GrandContigs) {
		push @Ordered, sprintf("%10.10s\t%s", $Leftmost{$GrandContigs[$i][0]}, $i);
	}
	@OkOverlap = ("r0144J01", ); %OkToOverlap = ();
	foreach $bac (@OkOverlap) {
		$OkToOverlap{$bac} = "Y";
	}
	@Ordered = sort {$a cmp $b} @Ordered;
	my $chunkCounter = 0;
	foreach $o (@Ordered) {
		my $assembled = ""; my $overhang = ""; $assemLength = 0;
		($dumb,$i) = split "\t", $o;
		print "<FONT COLOR=red><B>$features{$GrandContigs[$i][0]}</B></FONT>\n" if (exists $features{$GrandContigs[$i][0]});
		foreach $bac (@{$GrandContigs[$i]}) {
			my @Sequences = my @Names = ();
			my $numLines = &LOADFASTA("$dir/$bac", \@Sequences, \@Names);
			my $seq = $Sequences[0];
			my $length = length ($seq);
			my $copyStart = 0;
			if ($ButtJoin{$lastBac} = $bac) {
				$copyStart = 6;				# Overlap is only the restriction site
			} elsif (($assemLength = length($assembled)) > 0) {
				$overhang		= substr($assembled, $assemLength-$biteSize-1);
				$copyStart		= index($seq, $overhang) + length($overhang);
				my $leftOver	= substr($assembled, $assemLength-$copyStart);
				my $rightOver	= substr($seq, 0, $copyStart);
				$fail = 0;
				if ($leftOver ne $rightOver) {
					$fail = 1;
					print "<FONT COLOR=red>";
					unless (exists $OkToOverlap{$bac}) {
						print "\nOVERLAP FAILURE for BAC $bac\n</FONT>";
						print "$leftOver\n$rightOver\n";
						die;
					}
				}
			} else {
				$copyStart = 0;
			}
			$assembled .= substr($seq,$copyStart);
			printf ("%10.10s ", $bac);
			print "</FONT>" if ($fail);
			$lastBac = $bac;
		}
		print "\n";
		$chunkCounter++;
		my $blockSize = 80;
		my $nameTack = $chunkCounter; $nameTack = "0" . $nameTack if ($nameTack < 10);
		my $fName = sprintf("Chunk%s", $nameTack);
		open (FILE, ">$HOME/Data/ChrYChunks/$fName") or die "LOADFILES: Can't write to $fName: $!\n";
		my $pos = 0;
		$assemLength = length($assembled);
		print FILE ">ChrY Sequence $fName\n";
		while ($pos < $assemLength) {
			print FILE substr($assembled, $pos, $blockSize);
			print FILE "\n";
			$pos += $blockSize;
		}
		close FILE;
		system ("chmod 666 $HOME/Data/ChrYChunks/$fName");
	}
	print "</PRE>";
}

# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PRINTOPTIONS {
	print "<TABLE BORDER=1>";
	$colwidth = 250;
	print "<TR><TD ALIGN=center COLSPAN=2 BGCOLOR=GOLD><B><FONT SIZE=+2>Super-Mega-Mega Dot Plot</FONT></B><BR>";
	
	print "<TR><TD VALIGN=TOP WIDTH=$colwidth>";
	print $query->start_multipart_form(
								-action=>"Dotplot.cgi",
								-method=>POST,
								-name=>"MainForm");
	print "<FONT COLOR=darkgreen SIZE=+1><B>X-Axis Sequence</B></FONT><BR>\n";
	print $query->filefield(	-name=>'uploaded_file',
								-default=>$filename,
								-title=>"A text file containing a single DNA sequence in FastA format. Use the Browse... button to select the file, do not type it in directly.",
								-size=>20),"<BR>\n";
	print "<FONT COLOR=orange SIZE=+1><B>Y-Axis Sequence</B></FONT><BR>\n";
	print $query->filefield(	-name=>'other_file',
								-title=>"A text file containing supplemental information. Currently only Genescan output files can be interpreted.",
								-default=>$othername,
								-size=>20),"<BR>\n";
	print $query->submit(		-value=>'Analyze the indicated files',
								-title=>"Run the program with the settings and files indicated above.",
								),"\n";
	print "<TD VALIGN=TOP WIDTH=$colwidth>\n";
	print "Basepairs per pixel edge ";
	print $query->textfield(	-name=>'pixSize',
								-title=>"The number of bases that each pixel will span.",
								-default=>$pixSize,
								-size=>4),"<BR>\n";
	print $query->endform;
	print "</TABLE>";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Netscape does not clean up hex representations of special characters (ie leaves things like ":" as "%3A")
sub TIDYBROWSER {
	undef (@tidy);
	@tidy = split "%", $_[0];	
	$cleaned = $tidy[0];
	for $k (1..@tidy) {
		$n = substr ($tidy[$k],0,2);
		$cleaned .= chr(hex($n)) . substr($tidy[$k],2);
	}
	return $cleaned;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
#	[0] = the sequence to be analyzed	[1] = the number of bases in the loop of a hairpin
sub COLORSEQ {
	$a = $_[0];
	$l = (length($a) - $_[1]) / 2;	# number of bases on either side of loop
	$out = ""; $c = 0;
	while ($a ne "") {
		$_ = chop $a; $c++;
		if ($_[1] && ($c > $l) && ($c <= $l + $_[1])) {
			$out = "<FONT COLOR=purple SIZE=-2>$_</FONT>" . $out;
		}elsif (/a/i) {
			$out = "<FONT COLOR=green>A</FONT>" . $out;
		} elsif (/c/i) {	
			$out = "<FONT COLOR=blue>C</FONT>" . $out;
		} elsif (/t/i) {	
			$out = "<FONT COLOR=red>T</FONT>" . $out;
		} else {	
			$out = $_ . $out;
		}
		
	}
	return $out;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub SETCOLOR {
	$white		= $gph->colorAllocate(255,255,255);
	$black		= $gph->colorAllocate(0,0,0);
	$gph->transparent($white);
	$gph->interlaced('true');
	$gray50		= $gph->colorAllocate(128,128,128);
	$gray75		= $gph->colorAllocate(192,192,192);
	$fog		= $gph->colorAllocate(92,186,201);
	$watercress	= $gph->colorAllocate(187,212,26);
	$fern		= $gph->colorAllocate(92,163,31);
	$marigold	= $gph->colorAllocate(255,204,23);
	$orange		= $gph->colorAllocate(255,145,23);
	$melon		= $gph->colorAllocate(255,117,23);
	$salmon		= $gph->colorAllocate(234,97,87);
	$carnation	= $gph->colorAllocate(217,69,107);
	$apple		= $gph->colorAllocate(234,31,28);
	$orchid		= $gph->colorAllocate(186,87,194);
	$tan		= $gph->colorAllocate(191,133,74);
	$dirt		= $gph->colorAllocate(140,97,54);
	$lemmon		= $gph->colorAllocate(255,235,23);
	
	undef (@spect); # This is the color spectrum used by TNG RH
	$spect[0]	= $white;
	$spect[1]	= $gph->colorAllocate(214,214,214);
	$spect[2]	= $gph->colorAllocate(189,189,189);
	$spect[3]	= $gph->colorAllocate(156,156,156);
	$spect[4]	= $gph->colorAllocate(0,255,206);
	$spect[5]	= $gph->colorAllocate(0,255,90);
	$spect[6]	= $gph->colorAllocate(33,255,0);
	$spect[7]	= $gph->colorAllocate(148,255,0);
	$spect[8]	= $gph->colorAllocate(255,239,0);
	$spect[9]	= $gph->colorAllocate(255,123,0);
	$spect[10]	= $gph->colorAllocate(255,0,0);

	
	$red		= $gph->colorAllocate(255,0,0);
	$green		= $gph->colorAllocate(0,128,0);
	$yellow		= $gph->colorAllocate(255,255,0);
	$blue		= $gph->colorAllocate(0,0,255);
	$lightblue		= $gph->colorAllocate(128,128,255);
	for ($k=20; $k <= 100; $k += 5) {
		$j = int(255 * (100-$k) / 100);
		$blue{$k} = $gph->colorAllocate($j,$j,255);
	}
	$blue		= $gph->colorAllocate(0,0,255);
	$darkred	= $gph->colorAllocate(128,0,0);
	
###	@score =	(0,			30,			50,			100,		200,		300,			350,		400,		450,		500);
#	@ScoreCol	= ($white,	$gray75,	$gray50,	$fog,		$fern,		$watercress,	$marigold,	$melon,		$apple,		$orchid,	$green);
	@ScoreCol	= @spect;
	@RepeatCol	= ($white,	$blue{20},	$blue{30},	$blue{40},	$blue{50},	$blue{60},		$blue{70},	$blue{80},	$blue{90},	$blue{95},	$blue{100});
	if ($ColorType) {
		for $i (0..$#ScoreCol) { $ScoreCol[$i] = $gray50; }
		@TagCol = ($marigold, $lightblue, $melon, $red, $lemmon, $orchid, $fern, $watercress, $blue, $tan, $dirt);
		$TagCol[$TagCounter] = $dirt;
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub REVERSE {
	$d = $_[0]; $rev = "";
	while (length($d) > 0) {
		$n = chop $d;
		if ($rc{$n}) {
			$rev .= $rc{$n};
		} else {
			$rev .= $n;
		}
	}
	return $rev;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub MEGADOT {
	my $xRef	= shift;	# [0] Reference to string holding xAxis sequence
	my $yRef	= shift;	# [1] Reference to string holding yAxis sequence
	my $xName	= shift;	# [2] Name of sequence on x axis
	my $yName	= shift;	# [3] Name of sequence on y axis
	my $word	= shift;	# [4] Size in base pairs of words to be analyzed
	my $pixSize	= shift;	# [5] Number of base pairs that the edge of one pixel spans
	
	my $ratio = 1 / $pixSize;
	my $xLen = length($$xRef);
	my $yLen = length($$yRef);
	
	my $wT	= (gdTinyFont->width);			# Width in pixels of one character
	my $hT	= (gdTinyFont->height);			# Width in pixels of one character
	my $wM	= (gdMediumBoldFont->width);
	
	$width	= int($xLen/$pixSize);
	$height	= int($yLen/$pixSize);	
	$border	= $wT * 10;
	
	$gph = new GD::Image($width + 2 * $border, $height + 2 * $border);	
	$white		= $gph->colorAllocate(255,255,255);
	$gray95		= $gph->colorAllocate(242,242,242);
	$gray90		= $gph->colorAllocate(225,225,225);
	$gray75		= $gph->colorAllocate(192,192,192);
	$gray65		= $gph->colorAllocate(166,166,166);
	$gray50		= $gph->colorAllocate(128,128,128);
	$gray25		= $gph->colorAllocate(64,64,64);
	$black		= $gph->colorAllocate(0,0,0);
	$ltblue		= $gph->colorAllocate(100,100,255);
	$blue		= $gph->colorAllocate(0,0,255);
	$ltgreen	= $gph->colorAllocate(0,255,0);
	$green		= $gph->colorAllocate(0,128,0);
	$yellow		= $gph->colorAllocate(255,255,0);
	$orange		= $gph->colorAllocate(255,128,0);
	$red		= $gph->colorAllocate(255,0,0);
	$fog		= $gph->colorAllocate(92,186,201);
	$watercress	= $gph->colorAllocate(187,212,26);
	$fern		= $gph->colorAllocate(92,163,31);

	my $prob	= 4; for $i (2..$word) { $prob *= 4; }		# $prob is ($word)^4
	my $pixprob	= int(2*$pixSize*$pixSize / $prob) + 1;		# The number of hits expected at random in a given pixel (factor of 2 includes the reverse frame)
	@Type =		("Match Class",	"Random",	"2x Random",	"1/4 Diagonal", "1/2 Diagonal",	"3/4 Diagonal", "Diagonal", "Over Diagonal"	);	# Text Labels for scoring scheme
	@CutCol =	($black,		$gray95,	$gray90,		$fog,			$watercress,	$fern,			$orange,	$red			);	# Color Assignments for plotting
	@WebCol =	("black",		"gainsboro","darkgray",		"skyblue",		"greenyellow",	"forestgreen",	"orange",	"red"			);	# Color for use in HTML code
	@CutOff =	(0,				$pixprob,	$pixprob*2);																					# Hit cutoffs for that particular color
	$a = $CutOff[1]; $Hits = int (0.99 + $pixSize/2) - $a - 1;																				# The number of hits between randomness and the diagonal
	$Hits /= 3;																																# Adjust divisor (and above labels) for number of categories desired
	push @CutOff,											( int($Hits+0.5)+$a,int(2*$Hits+0.5)+$a,int(3*$Hits)+$a);						# Cutoffs for the intermediate hits
	push @CutOff,																							$pixSize + $a;					# NOTE: Diagonal, eg the match expected along the diagonal, starts at HALF the edge size. This is because if a diagonal-class match bisects a pixel halfway along the edge, the match will be HALF the hits received for a perfect bisection. Added padding of 2xpixprob is added at the end
	push @CutOff,																										$pixSize * 10;		# The upper limit of plotability - values outside this will not be plotted!
	
	$col[0] = $white;								# Color for no hit (should just be skipped anyway in plotting routine)
	$maxW = 0; $maxOverall = 0;
	for $i (1..$#Type) {
		$first = $CutOff[$i-1]+1;					# The previous cutoff value +1
		for $j ($first..$CutOff[$i]) {				# Cycle through all values for this type, assign the proper color to each
			$col[$j] = $CutCol[$i];
		}
		$maxW = length($Type[$i]) if (length($Type[$i]) > $maxW);
		$Range[$i] = " " . $first;
		$Range[$i] .= "-" . $CutOff[$i] if ($first != $CutOff[$i]);
		$maxOverall = length($Type[$i]) + length($Range[$i]) if (length($Type[$i]) + length($Range[$i]) > $maxOverall);
	}
	$Range[0] = " #Hits";
	
	$_ = $xName; s/>//;
	$x = $border + ($width - $wM*length($_))/2;
	$gph->string(gdMediumBoldFont,$x, 0, $_, $orange);				# Put name of Primary sequence on the x axis
	$_ = $yName; s/>//;
	$y = $border + ($height + $wM*length($_))/2;
	$gph->stringUp(gdMediumBoldFont,0, $y, $_, $orange);			# Put name of Secondary sequence on the y axis
	
	my $lineCounter = -1;
	$SquareKb = $xLen * $yLen / 1000000;											# The Area of the plot, in square kilobases
	
	my $maxTicks = $width / ($hT + 2);							# The maximum number of x-axis values that can legibily be displayed
	my $tickStep = int ($xLen / $maxTicks);					# The spacing that such ticks would have
	my @NiceTicks = (1,2,5,10);
	my $factor = 1; my $nt = 0; my $factorMod = "";
	while ($NiceTicks[$#NiceTicks] * $factor < $tickStep) { $factor *= 10; }	# Find the order of magnitude that encompases the tickstep
	while ($NiceTicks[$nt] * $factor < $tickStep) { $nt++;}						# Find the appropriate "nicetick" that suits tickstep
	$factorMod = " kb" if ($factor >= 1000);
	for ($k = 0; $k <= $xLen; $k += $NiceTicks[$nt] * $factor) {					# Add annotation to the top x axis. The  bottom axis will be annotated by other subroutines (eg FastA searching)
		$x = $k * $ratio;
		my $tickLab = $k;
		$tickLab = $k / 1000 . $factorMod if ($factorMod ne "");
		my $siz = $wT*length($tickLab);
		$gph->line($x+$border,$border-3,$x+$border,$border,$black);								# Draw X ticks, top
		$gph->line($x+$border,$border+$height,$x+$border,$border+$height+3,$black);				# Draw X ticks, bottom
		$gph->stringUp(gdTinyFont,$border+$x-($hT/2),$border-6,$tickLab,$black);				# Draw X labels, top
		$gph->stringUp(gdTinyFont,$border+$x-($hT/2),$border+$height+6+$siz,$tickLab,$black);				# Draw X labels, bottom
		$lineCounter++; next if ($lineCounter % 5);												# Draw only every 5th gridline
		$gph->line($border+$x,$border,$border+$x,$border+$width,	$yellow);					# Draw X Gridlines
	}
	$lineCounter = -1;
	for ($k = 0; $k <= $yLen; $k += $NiceTicks[$nt] * $factor) {								# Add annotation to both Y axes
		$x = $k * $ratio;
		my $tickLab = $k;
		$tickLab = $k / 1000 . $factorMod if ($factorMod ne "");
		my $siz = $wT*length($tickLab);
		$gph->line($border-3,$x+$border,$border,$x+$border,$black);								# Draw Y ticks, left side
		$gph->line($border+$width,$x+$border,$border+$width+3,$x+$border,$black);				# Draw Y ticks, right side
		$gph->string(gdTinyFont,$border-6-$siz, $x-($hT/2)+$border,$tickLab,$black);	# Draw Y labels, left side
		$gph->string(gdTinyFont,$border+$width+6, $x-($hT/2)+$border,$tickLab,$black);			# Draw Y labels, right side
		$lineCounter++; next if ($lineCounter % 5);												# Draw only every 5th gridline
		$gph->line($border+1,$x+$border,$border+$width,$x+$border,$yellow);						# Draw Y Gridlines
	}
	
	undef(%yHash);
	for ($y = 0; $y < $yLen - $pixSize - $word; $y += $pixSize) {
		undef(%{$yHash{$y}});
		for ($w = $y; $w < $y + $pixSize; $w++) {
			my $oligo = substr($$yRef, $w, $word);
			$yHash{$y}{$oligo}++;
			$yHash{$y}{&REVERSE($oligo)}++;
		}
	}
	for ($x = 0; $x < $xLen - $pixSize - $word; $x += $pixSize) {
		undef(%xHash);
		for ($w = $x; $w < $x + $pixSize; $w++) {
			my $oligo = substr($$xRef, $w, $word);
			$xHash{$oligo}++;
			$xHash{&REVERSE($oligo)}++;
		}
		for ($y = 0; $y < $yLen - $pixSize - $word; $y += $pixSize) {
			my $total = 0;
			foreach $n (keys %xHash) {
				$total += $xHash{$n} * $yHash{$y}{$n};
			}
			next if ($total < 1);
			$gph->setPixel($border+$x/$pixSize,$border+$y/$pixSize,$col[$total]);			# Set pixel to appropriate color
		}
	}

	$keyX = $border + $width - $wT * $maxOverall;
	$keyY = $border + 1;
	$gph->filledRectangle($keyX,$border,$width+$border, $border + ($hT+1) * ($#Type+2) + 1, $white);	# blank out image in area where key will be
	for $i (0..$#Type) {
		$a = sprintf ("%${maxW}s",  $Type[$i]); $a .= $Range[$i];
		$gph->string(gdTinyFont,$keyX, $keyY, $a, $CutCol[$i]);	# Draw Y labels
		$keyY += $hT+1;
	}
	$gph->string(gdTinyFont,$keyX, $keyY, "$pixSize bp/px, $word bp Word", $blue);	# Indicate pixel size and word size on key

	
	return $gph;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub DOTPLOT {
	$DotTime = time;
	undef(%WordHash); undef(%BackHash);
	$PixelStep = $width / $NumberOfBP;							# Pixels per base pair, used to determine pixel bin in later loops
	$halfDot = int(($dotWord-1)/2);
	for $i ($halfDot..($NumberOfBP - 1 - $halfDot)) {
		$a = substr ($QuerySeq, $i-$halfDot, $dotWord);
		push @{$WordHash{$a}}, $i;
	}
	
	$pixelEdge = int($NumberOfBP / $width) + 1;					# The number of bases spaning the edge of a pixel.
	$prob = 4; for $i (2..$dotWord) { $prob *= 4; }				# $prob is ($dotWord)^4
	$pixprob = int(2*$pixelEdge*$pixelEdge / $prob) + 1;		# The number of hits expected at random in a given pixel (factor of 2 includes the reverse frame)
	$top = $left;
	$height = $width;
	undef (@Matirx); $max = 0;									# Tally up the hits in each pixel
	$SimularityIndex = $RevSimIndex = 0;
	
	if ($DotBP) {													# A secondary sequence has been specified
		$height = $width * ($DotBP/$NumberOfBP);					# Respecify the height of the image to accomodate the size of the secondary sequence
		for $i ($halfDot..($DotBP - 1 - $halfDot)) {
			$a = substr ($DotSeq, $i-$halfDot, $dotWord);
			push @{$BackHash{$a}}, $i;
		}
		foreach $k ( keys %WordHash ) {								# Cycle through the identified words
			foreach $x ( @{$WordHash{$k}}) {						# Cycle along x-axis for each word
				$plotX = int ($PixelStep * $x);
				foreach $y ( @{$BackHash{$k}}) {					# Cycle along the y-axis for each word match in the secondary sequence
					$plotY = int ($PixelStep * $y);
					$Matrix[$plotX][$plotY]++;						# Increment the appropriate pixel
				}
				@Rev = @{$BackHash{&REVERSE($k)}};
				foreach $y ( @Rev ) {								# Cycle along the reverse complement of the secondary sequence
					$plotY = int ($PixelStep * $y);					
					$Matrix[$plotX][$plotY]++;						# Increment the appropriate pixel
				}
				$SimularityIndex	+= ($#{$WordHash{$k}}+1) * ($#{$BackHash{$k}}+1);
				$RevSimIndex		+= ($#{$WordHash{$k}}+1) * ($#Rev + 1);
			}
		}
	} else {														# The two cases (self-compare vs. two sequence compare) have been put in totally separate loops in the interest of speed
		foreach $k ( keys %WordHash ) {								# Cycle through the identified words
			@For = @{$WordHash{$k}};
			foreach $x ( @For ) {									# Cycle along x-axis for each word
				$plotX = int ($PixelStep * $x);
				foreach $y ( @For ) {								# Cycle along the y-axis for each word against itself
					$plotY = int ($PixelStep * $y);
					$Matrix[$plotX][$plotY]++;						# Increment the appropriate pixel
				}
				@Rev = @{$WordHash{&REVERSE($k)}};
				foreach $y ( @Rev ) {								# Cycle just the Y axis - the X axis will be cycled when the actual complement is reached
					$plotY = int ($PixelStep * $y);					# 
					$Matrix[$plotX][$plotY]++;						# Increment the appropriate pixel
				}
			}
			$SimularityIndex	+= ($#For+1) * ($#For+1);
			$RevSimIndex		+= ($#For+1) * ($#Rev+1);
		}
	}
	undef(%WordHash); undef(%BackHash);			# No longer need the hashes, free up the memory
	$DotTime = time - $DotTime;
	$DrawTime = time;
	
	$gph = new GD::Image($width + 2 *$left,$height + $top);	
	$white = $gph->colorAllocate(255,255,255);
	$gray95 = $gph->colorAllocate(242,242,242);
	$gray90 = $gph->colorAllocate(225,225,225);
	$gray75 = $gph->colorAllocate(192,192,192);
	$gray65 = $gph->colorAllocate(166,166,166);
	$gray50 = $gph->colorAllocate(128,128,128);
	$gray25 = $gph->colorAllocate(64,64,64);
	$black = $gph->colorAllocate(0,0,0);
	$ltblue = $gph->colorAllocate(100,100,255);
	$blue = $gph->colorAllocate(0,0,255);
	$ltgreen = $gph->colorAllocate(0,255,0);
	$green = $gph->colorAllocate(0,128,0);
	$yellow = $gph->colorAllocate(255,255,0);
	$orange = $gph->colorAllocate(255,128,0);
	$red = $gph->colorAllocate(255,0,0);
	$fog		= $gph->colorAllocate(92,186,201);
	$watercress	= $gph->colorAllocate(187,212,26);
	$fern		= $gph->colorAllocate(92,163,31);
	$gph->transparent($white);
	$gph->interlaced('true');
	$grays = 100;
	$DotSpectrum = 1;
	
	$wM = (gdMediumBoldFont->width);
	$_ = $QueryName; s/>//; $x = $left + ($width - $wM*length($_))/2;
	$gph->string(gdMediumBoldFont,$x, 0, $_, $orange);				# Put name of Query sequence on the x axis
	if ($DotBP) { $_ = $DotName; s/>//; }
	$y = $top + ($height + $wM*length($_))/2;
	$gph->stringUp(gdMediumBoldFont,0, $y, $_, $orange);			# Put name of secondary sequence on the y axis
	
	$maxTicks = $width / ($hT +2);	# The maximum number of x-axis values that can legibily be displayed
	$tickStep = int ($NumberOfBP / $maxTicks);				# The spacing that such ticks would have
	@NiceTicks = (1,2,5,10);
	$factor = 1; $nt = 0; $factorMod = "";
	while ($NiceTicks[$#NiceTicks] * $factor < $tickStep) { $factor *= 10; }	# Find the order of magnitude that encompases the tickstep
	while ($NiceTicks[$nt] * $factor < $tickStep) { $nt++;}						# Find the appropriate "nicetick" that suits tickstep
	$factorMod = " kb" if ($factor >= 1000);
	$HowMany = $NumberOfBP; my $lineCounter = -1;
	$SquareKb = $NumberOfBP * $NumberOfBP / 1000000;								# The Area of the plot, in square kilobases
	for ($k = 0; $k <= $HowMany; $k += $NiceTicks[$nt] * $factor) {					# Add annotation to the top x axis. The  bottom axis will be annotated by other subroutines (eg FastA searching)
		$x = $k * $ratio;
		my $tickLab = $k;
		$tickLab = $k / 1000 . $factorMod if ($factorMod ne "");
		$gph->line($x+$left,$top-3,$x+$left,$top,$black);							# Draw X ticks, top
		$gph->stringUp(gdTinyFont,$left+$x-($hT/2),$top-6,$tickLab,$black);			# Draw X labels, top
		$lineCounter++; next if ($lineCounter % 5);									# Draw only every 5th gridline
		$gph->line($left+$x,$top,$left+$x,$top+$width,	$yellow);					# Draw X Gridlines
	}
	if ($DotBP) {
		$HowMany = $DotBP;															# Adjust the vertical scale to the secondary sequence
		$SquareKb = $DotBP * $NumberOfBP / 1000000;
	}
	$RandomHits = int(0.5+($HowMany-$dotWord+1) * ($NumberOfBP-$dotWord+1)/$prob);	# The number of hits expected if the sequence was perfectly random
	$PerfectHits = $HowMany; $PerfectHits = $NumberOfBP if ($NumberOfBP<$HowMany);	# Number of hits expected for a perfect match (and no duplicate words in the match)
	$lineCounter = -1;
	for ($k = 0; $k <= $HowMany; $k += $NiceTicks[$nt] * $factor) {					# Add annotation to both Y axes
		$x = $k * $ratio;
		my $tickLab = $k;
		$tickLab = $k / 1000 . $factorMod if ($factorMod ne "");
		$gph->line($left-3,$x+$top,$left,$x+$top,$black);							# Draw Y ticks, left side
		$gph->line($left+$width,$x+$top,$left+$width+3,$x+$top,$black);				# Draw Y ticks, right side
		$gph->string(gdTinyFont,$left-6-$wT*length($k), $x-($hT/2)+$top,$tickLab,$black);	# Draw Y labels, left side
		$gph->string(gdTinyFont,$left+$width+6, $x-($hT/2)+$top,$tickLab,$black);			# Draw Y labels, right side
		$lineCounter++; next if ($lineCounter % 5);									# Draw only every 5th gridline
		$gph->line($left+1,$x+$top,$left+$width,$x+$top,$yellow);					# Draw Y Gridlines
	}
	
	@Type =		("Match Class",	"Random",	"2x Random",	"1/4 Diagonal", "1/2 Diagonal",	"3/4 Diagonal", "Diagonal", "Over Diagonal"	);	# Text Labels for scoring scheme
	@CutCol =	($black,		$gray95,	$gray90,		$fog,			$watercress,	$fern,			$orange,	$red			);	# Color Assignments for plotting
	@WebCol =	("black",		"gainsboro","darkgray",		"skyblue",		"greenyellow",	"forestgreen",	"orange",	"red"			);	# Color for use in HTML code
	@CutOff =	(0,				$pixprob,	$pixprob*2);																					# Hit cutoffs for that particular color
	$a = $CutOff[1]; $Hits = int (0.99 + $pixelEdge/2) - $a - 1;																		# The number of hits between randomness and the diagonal
	$Hits /= 3;																															# Adjust divisor (and above labels) for number of categories desired
	push @CutOff,											( int($Hits+0.5)+$a,int(2*$Hits+0.5)+$a,int(3*$Hits)+$a);								# Cutoffs for the intermediate hits
	push @CutOff,																							$pixelEdge + $a;				# NOTE: Diagonal, eg the match expected along the diagonal, starts at HALF the edge size. This is because if a diagonal-class match bisects a pixel halfway along the edge, the match will be HALF the hits received for a perfect bisection. Added padding of 2xpixprob is added at the end
	push @CutOff,																										$pixelEdge * 10;	# The upper limit of plotability - values outside this will not be plotted!
	
	$col[0] = $white;								# Color for no hit (should just be skipped anyway in plotting routine)
	$maxW = 0; $maxOverall = 0;
	for $i (1..$#Type) {
		$first = $CutOff[$i-1]+1;					# The previous cutoff value +1
		for $j ($first..$CutOff[$i]) {				# Cycle through all values for this type, assign the proper color to each
			$col[$j] = $CutCol[$i];
		}
		$maxW = length($Type[$i]) if (length($Type[$i]) > $maxW);
		$Range[$i] = " " . $first;
		$Range[$i] .= "-" . $CutOff[$i] if ($first != $CutOff[$i]);
		$maxOverall = length($Type[$i]) + length($Range[$i]) if (length($Type[$i]) + length($Range[$i]) > $maxOverall);
	}
	$Range[0] = " #Hits";
	$keyX = $left + $width - $wT * $maxOverall;
	$keyY = $top + 1;
	
	undef (@Hist );		# Histogram of hits
	for $x (0..$width) {
		for $y (0..$height) {
			next unless ($Matrix[$x][$y]);
			$Hist[ $Matrix[$x][$y] ]++;										# Increment the histogram array
			$gph->setPixel($x+$left,$y+$top,$col[$Matrix[$x][$y]]);			# Set pixel to appropriate color
		}
	}
	for $i (1..$#Hist) { $TotHist += $Hist[$i]; }		# Find the total number of pixels hit
	$Hist[0] = $width * $height - $TotHist++;			# Total number of pixels, minus those with hits, equals hitless pixels
	
	$gph->filledRectangle($keyX,$top,$width+$left, $top + ($hT+1) * ($#Type+2) + 1, $white);	# blank out image in area where key will be
	for $i (0..$#Type) {
		$a = sprintf ("%${maxW}s",  $Type[$i]); $a .= $Range[$i];
		$gph->string(gdTinyFont,$keyX, $keyY, $a, $CutCol[$i]);	# Draw Y labels
		$keyY += $hT+1;
	}
	$gph->string(gdTinyFont,$keyX, $keyY, "$pixelEdge bp/px, $dotWord bp Word", $blue);	# Indicate pixel size and word size on key
	
	$Dotname = "Dotgraph$$" . ".tmp";
	open (FILE, ">$DataLoc/$Dotname") or die "Failure to write to $gphname: $!\n";
	print FILE $gph->gif;
	close FILE;
	$DrawTime = time - $DrawTime;
}
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
sub DEFINECONSTANTS {
	@Environmentals = (	HTTP_FROM, REMOTE_USER, SERVER_NAME, REMOTE_HOST, HTTP_USER_AGENT,
						REMOTE_ADDR, HTTP_REFERER, SERVER_SOFTWARE, SERVER_PROTOCOL,
						AUTH_TYPE, REMOTE_IDENT);
	$HOME = "/usr/local/etc/httpd/cgi-bin/pagelab/charles";
	$Here = "."."$ENV{'REMOTE_ADDR'}";				# Returns the IP address of the *USER'S* Machine (to allow multiple users)
	
	$DataLoc =	"/tmp/ChasFasta" . $Here;			# Unique directory for the user in the tmp folder
	$WebLoc =	"/tmp/ChasFasta" . $Here;			# The base URL for the host computer
	unless (-e $DataLoc) {							# Skip this if the directory already exists
		system ("mkdir $DataLoc");
		print "<FONT COLOR=brick>Temporary folder $DataLoc created for data storage.</FONT><BR>\n";
		system("chmod 777 $DataLoc");
	}
	
	@RandNucs	= ("A","C","T","G",);
	$Allowed	= "AGCTUYRWSKMBDHVXN";
	$Reverse	= "TCGAARYWSMKVHDBXN";
	%Ambig 		= ( A => 'A',	M => 'AC',	S => 'CG',	V => 'ACG',	X => 'ACGT',
					C => 'C',	R => 'AG',	Y => 'CT',	H => 'ACT',	N => 'ACGT',
					T => 'T',	W => 'AT',	K => 'GT',	D => 'AGT',
					G => 'G',							B => 'CGT', 				);
	
	undef (%twoAmbig); undef (%RevAmbig);
	foreach my $a (keys %Ambig) {
		$RevAmbig{$Ambig{$a}} = $a;				# Hash that points in the reverse of Ambig
	}
	foreach my $c (keys %Ambig) {
		foreach my $d (keys %Ambig) {
			undef (@t); undef (%u);
			@s = split '', $Ambig{$c};			# @s will contain all the possible bases for this *pair* of ambiguities
			push @s, (split '', $Ambig{$d});	# Load up @s with the possibilities
			foreach my $n (@s) {
				push @t, $n if ($u{$n} eq "");	# Load components into @t unless they're already there
				$u{$n} = "Y";					# %u keeps track of what has gone into @t already
			}
			@t = sort {$a cmp $b} @t;			# sort @t alphabetically
			my $j = join '', @t;				# $j represents a short string of all the bases
			$twoAmbig{$c}{$d} = $RevAmbig{$j};
#			print "($c = $Ambig{$c}) + ($d = $Ambig{$d}) = $j = $twoAmbig{$c}{$d}<BR>";
		}
	}
#	%rc = (	"A" => "T", C => G, T => A, G => C, );	# Previous way of defining reverse complement
	%features = (	"1188O08" => "GAP 1",		"0071C04" => "TSPY Gap",
					"0344D02" => "TSPY Gap",	"0075F05" => "Centromere",
					"0509B06" => "GAP 2",		"0268K13" => "GAP 3",
					"0400I17" => "GAP 4",		"0242E13" => "Heterochromatic Boundary",
					"0057J19" => "Long Arm PA",	);
	%rc = ();
	for $i (0..length($Allowed)-1) {
		$rc{ substr( $Allowed , $i ,1) } = substr( $Reverse , $i ,1 );
	}
	$BlockSize = 80;	# Number of nucs to print in a FastA block
	$hT = (gdTinyFont->height);
	$wT = (gdTinyFont->width);
	$minExonSize = 30;
	$minIntonSize = 30;
	$BigGap = "";
	for $i (1..$minIntonSize) {
		$BigGap .= "-";
	}
	
	$query = new CGI;
	$filename =	$query->param('uploaded_file');
	$libname =	$query->param('uploaded_lib');
	$othername = $query->param('other_file');
	
	@{$Rank{'Zscore'}}	= (0,	30,	50,	100,	200,	300,	350,	400,	450,	500);
	@{$Rank{'SWscore'}}	= (0,	30,	50,	100,	200,	300,	350,	400,	450,	500);
	@{$Rank{'PerID'}}	= (0,	50,	60,	70,		80,		85,		90,		95,		98,		100);
#	@{$Rank{'expect'}}	= (2,	0,	-2,	-4,		-8,		-12,	-14,	-18,	-25,	-30);
	
	@dataTypes = (	"Name",		"initn",	"init1",	"opt",			"Zscore",	"expect",
					"SWscore",	"PerID",	"Overlap",	"MatchStart",	"MatchEnd",	"LibStart",
					"LibEnd",	"Files",	);
	
	# The following routine retrieves values from the CGI paramater lists
	# It also sets defaults if no value is specified.				
	@CGIParams = (	"dotWord=15",		"pixSize=1000",);

	foreach $i (@CGIParams) {	
		($var, $val) = split "=", $i;
		$$var = $query->param($var) if ($query->param($var) ne "");
		$$var = $val if ($$var eq "");
	}
	
	#	$dotWord++ if (($dotWord/2) == int($dotWord/2));
	#	MinPalSize	=	Maximum reach to consider
	#	MaxPalSize	= 	Minimum reach to consider 0 = 2bp palindromes
	#	MinSimple	=	Minimum size to be considered a simple repeat
	#	SimpStep	=	Steps to take while looking for simple repeats (suggest 1 or 2)
	#	SimpSize	=	Maximum unit size of simple repeat to consider
	#	MinLoop		=	Minimum number of nucleotides to consider in the loop of a hairpin
	#	MaxLoop		=	Maximum number of nucleotides to consider in the loop of a hairpin
	#	MinHairpin	=	Minimum number of bases needed to be considered a hairpin
	#	BinSize		=	Size of window to analyze palindrome frequency
	#	ChunkSize	=	Maximum fragment size to pass to FastA
	$ChunkSize	= 20000 if ($ChunkSize > 20000);	# This is the maximum currently allowed by FastA
	$MinReach = $MinPalSize/2 - 1;
	$MaxReach = $MaxPalSize/2 - 1;
	
	# Minimum number of repeats to be considered valid
	@RepeatNames =	("Null", "Simple", "Dinucleotide", "Trinucleotide", "Tetranucleotide", "Pentanucleotide", "Hexanucleotide",	"Heptanucleotide");
	@SimpLimit =	(10000,		8,			5,				4,					3,				3,					3,					3);
	@F = (); 
	$citation = "FASTA citation: W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Prints a mini-plot to-scale with the x-axis ($width)
sub AUXGRAPH {
	my $maxRef		= shift;	# [0] - reference to data set A (maximum values)
	my $minRef		= shift;	# [1] - reference to data set B (minimum values, can be left out)
	my $plotHeight	= shift;	# [2] - height over which data should be plotted
	my $yLabel		= shift;	# [3] - Y axis label
	my $colRef		= shift;	# [4] - color to use for plotting
	my $tickStep	= shift;	# [5] - Minimum axis tick step
	
	my $Ctop= $hT; my $Cbottom = 3; my $height = $Ctop + $Cbottom + $plotHeight;
	my $gph	= new GD::Image($width + $left + $right, $height);
	$white	= $gph->colorAllocate(255,255,255);
	$blue	= $gph->colorAllocate(0,0,255);
	$green	= $gph->colorAllocate(0,255,0);
	$orange	= $gph->colorAllocate(255,128,0);
	$red	= $gph->colorAllocate(255,0,0);
	$black	= $gph->colorAllocate(0,0,0);
	$ltred	= $gph->colorAllocate(255,128,128);
	$gray75	= $gph->colorAllocate(192,192,192);
	$gph->transparent($white);
	$gph->interlaced('true');
			
	my $max = -10000000; my $min = 10000000;
	for $ag (0..$#{$maxRef}) {
		next if (${$maxRef}[$ag] eq "");
		$max = ${$maxRef}[$ag] if (${$maxRef}[$ag] > $max);
		$min = ${$maxRef}[$ag] if (${$maxRef}[$ag] < $min);
		if ($minRef != 0) {
			$max = ${$minRef}[$ag] if (${$minRef}[$ag] > $max);
			$min = ${$minRef}[$ag] if (${$minRef}[$ag] < $min);
		}
	}
	$min = int ($min / $tickStep) * $tickStep;
	$max = (1 + int (($max - 0.0000001) / $tickStep)) * $tickStep;

	$HalfHeight = ($height + $wT * length($yLabel) ) / 2;
	$gph->stringUp(gdTinyFont,$width+$left,$HalfHeight,$yLabel,${$colRef});
	$NumNiceTicks = 1; $k = $min;
	my $Vratio = $plotHeight / ($max - $min); 		# vertical plotting ratio
	while ( $NumNiceTicks * $tickStep * $Vratio  < $hT + 2) { $NumNiceTicks++; }
	while ($k <= $max) {
		$y =	$plotHeight + $Ctop - int(0.5 + ($k-$min) * $Vratio);
		$gph->line($left,$y,$width+$left,$y,$gray75);
		$gph->line($left-3,$y,$left,$y,$black);
		$d = " " . $k;
		$gph->string(gdTinyFont,$left-5 - $wT * length($d),$y-($hT/2),$d,${$colRef});
		$k += $tickStep * $NumNiceTicks;
	}
	for ($k = 0; $k <= $NumberOfBP; $k += $NiceTicks[$nt] * $factor) {	# Add X-grid
		$x = $k * $ratio + $left;
		$gph->line($x,$y,$x,$plotHeight + $Ctop,$gray75);
	}
	
#	$tag = "Window Size: "; $gph->string(gdTinyFont,$left,0,$tag,$green);
	my $xEx = $left + 10;
	my $xMin = int(($halfBin+1)*$ratio);
	my $xMax = int(($NumberOfBP-$halfBin)*$ratio);
	$gph->filledRectangle($xEx,1, $xEx + ($halfBin * 2 + 1) * $ratio, 4, $green);
	my $y1 = my $y2 = $plotHeight + $Ctop - int(0.5 + (${$minRef}[0]-$min) * $Vratio);
	for $xp ($xMin..$xMax) {
		if ($minRef != 0) {
			$y2 = $plotHeight + $Ctop - int(0.5 + (${$minRef}[$xp]-$min) * $Vratio);	# Min-Max lines
		} else {
			$y2 = $y1																	# Use the last plotted value
		}
		$y1 = $plotHeight + $Ctop - int(0.5 + (${$maxRef}[$xp]-$min) * $Vratio);
		next if ($y2 eq "");
		$gph->line($xp + $left,$y1,$xp + $left,$y2,${$colRef});
	}
	$_ = $yLabel; s/ //g; s/\///g; $_ = $_ . $$ . ".tmp";
	open (FILE, ">$DataLoc/$_") or die "Failure to write to $_: $!\n";
	print FILE $gph->gif;
	close FILE;
	print qq(<img src="$WebLoc/$_"><BR>\n);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
