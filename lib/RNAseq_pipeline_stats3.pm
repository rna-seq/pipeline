package RNAseq_pipeline_stats3;

# Export subroutines to caller namespace
use Exporter;
@ISA=('Exporter');
@EXPORT_OK=('get_Average','get_Stdev','get_Min_and_Max','get_Shannon_H',
	    'log10','get_Median','get_p_value_dist','get_factorial',
	    'get_p_value_hyp','get_t_significance','get_z_score',
	    'get_EffectSize');

use strict;
use warnings;

# Get the log base 10 of a number
sub log10 {
    my $n = shift;
    return log($n)/log(10);
}

# Get the effect size This is the difference of the averages divided by the
# Stdev in the control
sub get_EffectSize {
    my $control=shift;
    my $problem=shift;
    my $cstdev=shift;
    my $decimals=shift;
    $decimals=_set_decimals($decimals);

    if ($cstdev != 0) {
	return(sprintf("%.${decimals}f",abs($problem - $control) / $cstdev));
    } else {
	return('NAN');
    }
}

# Get the median of a string of values
sub get_Median {
    my $values=shift;
    my @sorted=sort {$a <=> $b} @{$values};
    my $size=@sorted;
    my $median;

    if ($size % 2) {
	$median=$sorted[int($size/2)];
    } else {
	$median=($sorted[($size / 2) - 1] + $sorted[$size / 2]) / 2;
    }
    return($median);
}

# Non recursive factorial function
sub get_factorial {
    my $number=shift;
    my $factorial=1;

    # account for nonintegers by rounding to the closest integer
    $number=int($number + 0.5);

    if ($number < 0) {
	die "Invalid number factorial of $number\n";
    } if ($number >= 0) {
	for (my $i=$number; $i > 1; $i--) {
	    $factorial*=$i;
	}
    }
    return($factorial);
}

# Get the z_score
sub get_z_score {
    my $obs=shift;
    my $avg=shift;
    my $stdev=shift;
    my $z_score=0;
    if ($stdev > 0) {
	$z_score=sprintf("%.6f",($obs - $avg)) / $stdev;
    }
    return($z_score);
}


# Determine if a t score is signeficant to 0.05 to 0.001 level in a t
# distribution with more than 1000 degrees of freedom in a two tailed
# distribution
sub get_t_significance {
    my $obs=shift;
    my $avg=shift;
    my $stdev=shift;
    my $score=abs(get_z_score($obs,
			      $avg,
			      $stdev));
     my $sig=1;
    
    if ($score > 3.29) {
	$sig=0.001;
    } elsif ($score > 2.80) {
	$sig=0.005;
    } elsif ($score > 2.58) {
	$sig=0.01;
    } elsif ($score > 1.96) {
	$sig=0.05;
    }
    return($sig);
}

# Get the binomial coefficient, which is the (x over y)
sub get_bin_coeff {
    my $x=shift;
    my $y=shift;
    my $diff=$x - $y;
    my $bin_coeff=(get_factorial($x)/
		   (get_factorial($y) * get_factorial($diff)));
    return($bin_coeff);
}
#### P Values
# These are several methods to obtain the p_value of an observed frequency
# Use the hypergeometric distribution formula to determine the probability of
# an expected observation
sub get_p_value_hyp {
    my $obs_successes=shift;
    my $pop_successes=shift;
    my $sample_size=shift;
    my $pop_size=shift;
    my $failures=$sample_size - $obs_successes;;
    my $pop_failures=$pop_size - $pop_successes;
    my $p_value=sprintf("%.6f",((get_bin_coeff($pop_successes,$obs_successes) *
				 get_bin_coeff($pop_failures,$failures)) /
				get_bin_coeff($pop_size,$sample_size)));
    return($p_value);
}
# Get p_value from the distribution directly (2 tailed or 1 tailed, one tailed
# by default). the minimum p_value it will return is 1e-10
# the reason for this is to allow for the p-value corrections for multiple
# testing
sub get_p_value_dist {
    my $values=shift;
    my $observed=shift;
    my $two_tailed=shift;
    unless($two_tailed) {
	$two_tailed=1;
    }
    my @sorted=sort {$a <=> $b} @{$values};
    my $size=@sorted;
    my $dist_from_low=0;
    my $dist_from_high=0;
    my $equal=0;
    my $median=get_Median($values);
    my $p_value;

    # Get the distance from the ends of the distribution
    for (my $i=0;$i<@sorted;$i++) {
	# from the lower end
	if ($sorted[$i] > $observed) {
	    $dist_from_low++;
	} elsif ($sorted[$i] == $observed) {
	    # If equal
	    $equal++;
	} elsif ($sorted[$i] < $observed) {
	    $dist_from_high++;
	} else {
	    warn "Problem\n";
	}
    }

    my $p_low=$dist_from_low / $size;
    my $p_high=$dist_from_high / $size;
    my $p_equal=$equal / $size;

#    print join("\t",
#	       $observed,
#	       $dist_from_low,
#	       $p_low,
#	       $dist_from_high,
#	       $p_high,
#	       'dist'),"\n";

    # If the test is two tailed we have to double the p_value as we don't know
    # in what direction the variation is.
    if ($two_tailed==2) {
	$p_low*=2;
	$p_high*=2;
	$p_equal*=2;
    }

    # Set a minimum, in case the number is 0, so we can use the p.adjust
    # from R
    if ($p_low == 0) {
	$p_low=1/$size;
    }
    if ($p_high== 0) {
	$p_high=1/$size;
    }

    if ($p_low <= $p_high) {
	$p_value=sprintf("%5g",$p_low + $p_equal);
    } else {
	$p_value=sprintf("%5g",$p_high + $p_equal);
    }
    
    return($p_value);
}

# Get the Shannon entropy of a string
sub get_Shannon_H {
    my $string=shift;
    my $decimals=shift;
    $decimals=_set_decimals($decimals);
    my %count;
    my $total;
    my $entropy=0;

    $string=uc($string);
    foreach my $char (split(//,$string)) {
	$count{$char}++;
	$total++;
    }

    foreach my $char (keys %count) {
	my $prob= $count{$char}/$total;
	$entropy+=$prob * log($prob);
    }

    return(sprintf("%.${decimals}f",-$entropy/log(2)));
}

# Get the average of an array of numbers
sub get_Average {
    my $values=shift;
    my $decimals=shift;
    $decimals=_set_decimals($decimals);
    my $count=@{$values};
    my $total=0;

    foreach my $value (@{$values}) {
	$total+=$value;
    }

    return(sprintf("%.${decimals}f",$total/$count));
}

sub get_Stdev {
    my $values=shift;
    my $decimals=shift;
    $decimals=_set_decimals($decimals);
    my $stdev;
    my $count=@{$values};
    my $mean=get_Average($values);
    foreach my $value (@{$values}) {
	$stdev+=( ($value - $mean) ** 2);
    }
    if ($count < 2) {
	warn "Unable to calculate Stdev with $count elements returning -1\n";
	return(-1);
    } else {
	return(sprintf("%.${decimals}f",sqrt($stdev/($count - 1))));
    }
}

sub get_Min_and_Max {
    my $values=shift;
    my @sorted=sort {$a <=> $b} @{$values};
    return($sorted[0],$sorted[-1]);
}

sub _set_decimals {
    my $decimals=shift;
    # Set the number of decimal places unless it is specified in the sub call
    unless ($decimals) {
	$decimals=2;
    }
    return($decimals);
}

1;
