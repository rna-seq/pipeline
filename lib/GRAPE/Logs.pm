package GRAPE::Logs;

use strict;
use warnings;

##################################################
## the object constructor (simplistic version)  ##
##################################################
sub new {
    my $class =shift;
    my $filename =shift;
    my $debug=shift;
	
    my $self  = {};

    $self->{debug} = $debug;
    $self->{logname}   = $filename;
    $self->{logfh}    = _get_fh($filename,$debug);
    $self->{start}  = time();
    $self->{end} = undef;
    $self->{duration}= 0;
    bless($self, $class);

    $self->printlog("Start: ",$self->{start});

    return $self;
}

sub printlog {
    my $self=shift;
    my $string=shift;

    print {$self->{logfh}} $string,"\n";

    return $self;
}

sub DESTROY {
    my $self=shift;
    $self->{end}=time();
    $self->{duration}=$self->{end} - $self->{start};
    printf {$self->{logfh}} ("\n\nTotal running time: %02d:%02d:%02d\n\n", 
	   int($self->{duration} / 3600),
	   int(($self->{duration} % 3600) / 60),
	   int($self->{duration} % 60));
    close($self->{logfh});
}

sub _get_fh {
    my $filename=shift;
    my $debug=shift;

    chomp($filename);
    my $openstring;
    my $fh;
    
    $openstring=">".$filename;

    if ($debug) {
	$fh=*STDERR;
    } else {
	open($fh,$openstring) ||
	    die "Unable to open $filename: $!,$?\n";
    }

    return($fh);
}

1;
