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
<html><TITLE>Staffa Process Killer</TITLE>
<BODY BGCOLOR=tomato  text="#000000" link="#118911">
EOF

require "ParseQuery.pl";
&PARSE;

if ($toKill) {
	$er = `cmd 2>&1`;
	$attempt = `kill -9 $toKill`;
	print "<H2>Attempt made to kill process $toKill.</H2>\n";
	print "<H2>Result: $attempt</H2>\n<H2>Status: $?</H2>\n<H2>Error: $er</H2>\n";
}

$output =	`ps -elf | grep perl`;
$output .=	`ps -elf | grep $toSearch` if ($toSearch);
$header = " F S      UID   PID  PPID  C PRI NI     ADDR     SZ    WCHAN    STIME TTY      TIME CMD\n";

print "<center><H1> The following processes were found:</H1></CENTER>\n";

print "<PRE><FONT COLOR=blue>$header</FONT>";
@temp = split "\n", $output;
foreach $_ (@temp) {
	next if (/grep/);
	print "$_\n";
}
print "<FONT COLOR=red>Keep in mind that at least one of the scripts above is *this* one that just executed.</FONT>\n";
print "</PRE>";

print <<EOF;
<FORM ACTION="ProcessKiller.cgi">
Also display processes that match:
 <INPUT TYPE=text NAME="toSearch" VALUE="$toSearch" SIZE=30><BR>
Kill this process: 
 <INPUT TYPE=text NAME="toKill" SIZE=6> (<B>KNOW WHAT YOU ARE DOING!</B><BR>
<INPUT TYPE=hidden NAME="RandNum" VALUE="$$">
</FORM>
EOF



