# Subroutine that loads fasta file, being on the lookout for junk
# from the resource forks of macintosh files
# Staffa path for this file is:
# /usr/local/etc/httpd/cgi-bin/pagelab/charles/

$| = 1;

sub LOADFASTA {
	my $fileTarget	= shift;	# [0] Filename
	my $seqRef		= shift;	# [1] Ref to array holding sequence information
	my $nameRef		= shift;	# [2] Ref to array holding names of the sequences
	
	my @lines = ();
	if ($fileTarget =~ /\//) {
		open (FILE, "$fileTarget") or die "LOADFASTA: Can't read $fileTarget: $!\n";
		while (<FILE>) {
			chomp;
			push @lines, (split /[\r\n]/, uc($_));	# Upper case all sequence
		}
		close FILE;
	} else {
		while (<$fileTarget>) {
			chomp;
			push @lines, (split /[\r\n]/, uc($_));	# Upper case all sequence
		}
	}
	$firstSeq	= 0;		# First line containing ONLY allowed sequence characters
	$lastSeq	= $#lines;	# Last line etc
	
	while ($firstSeq <= $lastSeq	&& &CHECKSEQ(\$lines[$firstSeq]))	{$firstSeq++;}
#	print "-> $firstSeq, $lastSeq<BR>\n";
	if ($firstSeq == 0) {
		&FASTAERR($fileTarget, "Could not find a '>' before the beginning of the sequence. Fasta files should be of the format<BR>>Sequence name<BR>(lines of sequence here)<BR>>Sequence name<BR>(lines of sequence here)<BR>etc.");
		return;
	}
	while ($firstSeq <= $lastSeq+1	&& &CHECKSEQ(\$lines[$lastSeq]))	{$lastSeq--;}
	
	if ($firstSeq > $lastSeq) {
		&FASTAERR($fileTarget, "Could not find any lines that contain pure sequence. Are you sure this is a Fasta file?");
		return;
	}
	$firstSeq--;
	$lastTag = rindex($lines[$firstSeq],">");			# Find the very last incidence of ">" before the first line of sequence
	$lines[$firstSeq] = substr($lines[$firstSeq], $lastTag);
	my $counter = -1;
	for $i ($firstSeq..$lastSeq) {
		my $theLine = $lines[$i];
		if (substr($theLine,0,1) eq ">") {				# Tag line indicating new sequence
			$counter++;
			${$nameRef}[$counter] = substr($theLine,1);
			next;
		}
		${$seqRef}[$counter] .= $theLine;
	}
	undef (@lines);
	return ($lastSeq - $firstSeq + 1);
}

sub CHECKSEQ {
	my $data = ${$_[0]};		# [0] is ref to string of sequence
#	print "$firstSeq = ..$data..<BR>\n";
	$data =~ s/[AGCTUYRWSKMBDHVXN]//g;
	if ($data eq "") {
		return 0;
	} else {
		return length($data);
	}
}

sub FASTAERR {
	print "<FONT COLOR=red>Error encountered while reading file $_[0]:</FONT><BR>\n";
	print "<FONT COLOR=brown>$_[1]</FONT><BR>\n";
	return;
print <<EOF;
<P>This annoyance generally occurs because some bizzare characters got included during a copy-and-paste operation. 
Certain programs written by Big Brother Bill will inapropriately load the file's resource fork as well as the data fork. 
This has the result of confusing the parsing algorithm. 
The best solution so far seems to be to <FONT COLOR=RED><B>use Netscape</B></FONT>. 
If you must still use Explorer, then try:<BR><BR>
<OL>
<LI>Copy the sequence only (not the name tag) out of your original FastA file.<BR>
<LI> Paste it into EditSeq.<BR>
<LI> Select all (command-A) and copy it.<BR>
<LI> Go to BBEdit, open a new document, and paste the sequence in.<BR>
<LI> Add a new FastA name tag (e.g. >MySequence)to the top, and save the file.<BR>
<LI> Try again with the new file.<BR>
</OL>
EOF
}
