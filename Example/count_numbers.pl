use strict;
use warnings;


use Getopt::Std;

my %OPTS;
getopts('i:o:nsm',\%OPTS);

my $file = $OPTS{"i"};
my $NUM_FILE;

open (NUM_FILE, "< $file") or die "Cannot open $file";

# Initialize a hash to store the counts of each number
my %number_counts;

# Read each line from the file
while (<NUM_FILE>) {
    chomp;

    # Split the line into individual numbers using whitespace as the delimiter
    my @numbers = split /\s+/, $_;

    # Count the occurrences of each number
    foreach my $number (@numbers) {
        $number_counts{$number}++;
    }
}

# Close the file
close $file;

# Print the counts
foreach my $number (sort keys %number_counts) {
    my $count = $number_counts{$number};
    print "Bottleneck $number appears $count times\n";
}
