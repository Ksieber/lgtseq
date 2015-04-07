#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/

=head1 NAME

lgtseq_download.pl

=head1 SYNOPSIS

Download a file from the TCGA, EGA, or SRA for lgtseq. 
This script can setup the pipeline to then launch lgtseq_prelim-filter.pl and/or lgtseq_analysis.pl

=head1 EXAMPLE

lgtseq_download.pl --input=84a23558-26b3-42f4-a9f3-d83d118f2327 --output_dir=/dir/for/output/ --qsub=1 --threads=4

=head1 AUTHOR - Karsten Sieber

e-mail: karsten.sieber@gmail.com

=cut

my $LGTSEQ_DL = '1.00';

use warnings;
no warnings 'uninitialized';
use strict;
use Carp;
$Carp::MaxArgLen = 0;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use print_call;
use lib ( "/local/projects-t3/HLGT/scripts/lgtseek/lib/", "/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/" );
use LGTSeek;
use File::Basename;
use setup_input;
use run_cmd;
use mk_dir;
use print_call;

my %options;
my $results = GetOptions(
    \%options,            'input|i=s',              'input_list|I=s',           'output_dir|o=s',          'threads|t=i',       'rate_limit=s',
    'cghub_key=s',        'Qsub|q=i',               'Qsub_iterate|Q=i',         'qsub_pause=i',            'sub_mail=s',        'sub_name=s',
    'verbose=i',          'sort_mem=s',             'sub_mem=s',                'tcga_dirs=i',             'subdirs=i',         'sub_mail=s',
    'project=s',          'launch_prelim_filter=i', 'prelim_dir=s',             'prelim_threads=i',        'prelim_dir=s',      'prelim_threads=i',
    'prelim_sort_mem=s',  'prelim_sub_mem=s',       'prelim_name_sort_input=i', 'prelim_delete_input=i',   'launch_analysis=i', 'analysis_threads=i',
    'analysis_sub_mem=s', 'analysis_dir=s',         'analysis_sort_mem',        'analysis_delete_input=i', 'help|h',            'help_full|?'
) or die "Error: Unrecognized command line option. Please try again.\n";

&help_full if ( $options{help_full} );    ## At the end of the script
&help if ( $options{help} );

# Check we have the correct inputs
if ( !$options{input} && !$options{input_list} ) { die "Must use --input\n"; }
if ( !$options{output_dir} ) { die "Must give use --output_dir=\"/dir/for/output/\"\n"; }

# Setup a few global variables
my $lgtseq              = LGTSeek->new2( \%options );
my $original_output_dir = $lgtseq->{output_dir};
my $inputs              = ( !$options{Qsub_iterate} ) ? setup_input( \%options ) : [ $options{input_list} ];

print_call( \%options, "LGTSEQ_DL_VERSION=$LGTSEQ_DL\tLGTSeek.pm_VERSION=$LGTSeek::VERSION" );

foreach my $input (@$inputs) {
    ## Setup output for each input
    $lgtseq->{tcga_dirs} = 1 if ( $input =~ /^\w{8}\-\w{4}\-\w{4}\-\w{4}\-\w{12}$/ or $input =~ /^\d{5,6}.*(\.f\w{0,3}q.gz)?$/ );
    my ( $fn, $path, $suf ) = fileparse( $input, ( @{ $lgtseq->{bam_suffix_list} }, @{ $lgtseq->{fastq_suffix_list} } ) );
    my $subdir = ( $suf =~ /\w+/ ) ? $fn : undef;
    my @split_path = split( /\//, $path );
    my $tcga_dir = ( $suf =~ /\.bam$|\.f\w{0,3}q$/ ) ? $split_path[-1] : $fn;
    $lgtseq->{output_dir} = $original_output_dir;
    if ( $lgtseq->{tcga_dirs} == 1 and !$options{Qsub_iterate} ) {
        $lgtseq->{output_dir} = $lgtseq->{output_dir} . "/$tcga_dir\/" unless ( $lgtseq->{output_dir} =~ /$tcga_dir\/*$/ );
    }
    if ( $lgtseq->{subdirs} == 1 and !$options{Qsub_iterate} and defined $subdir ) {
        $lgtseq->{output_dir} = $lgtseq->{output_dir} . "/$subdir\/" unless ( $lgtseq->{output_dir} =~ /$subdir\/*$/ );
    }
    $lgtseq->{output_dir} =~ s/\/{2,}/\//g;
    $options{output_dir} = $lgtseq->{output_dir};
    $lgtseq->_run_cmd("mkdir -p -m u=rwx,g=rwx,o= $lgtseq->{output_dir}");

    # Qsub this script foreach input and any of the options passed
    if ( $lgtseq->{Qsub} == 1 or $options{Qsub_iterate} == 1 ) {
        ## If we are in the orignal call, change input from list to a single file
        if ( $options{input_list} and !$options{Qsub_iterate} ) { $options{input} = $input; }
        ## Check $sub_mem is enough for sorting
        if ( $lgtseq->{name_sort_input} == 1 ) {
            my $original_sub_mem;
            my $original_sort_mem;
            if ( $lgtseq->{sub_mem} =~ /^(\d+)[K|M|G]$/ )  { $original_sub_mem  = $1; }
            if ( $lgtseq->{sort_mem} =~ /^(\d+)[K|M|G]$/ ) { $original_sort_mem = $1; }
            if ( $original_sub_mem < ( $original_sort_mem * $options{threads} ) ) {
                $options{sub_mem} = ( ceil( ( $original_sort_mem * $options{threads} ) * 1.1 ) ) + 1 . "G";
            }
        }
        ## Build qsub command
        my $cmd = "$^X $0";
        foreach my $key ( keys %options ) {
            next if ( $key eq 'Qsub' );
            next if ( $key eq 'Qsub_iterate' );
            next if ( $key eq 'qsub_pause' );
            next if ( $key eq 'sub_name' );
            next if ( $key eq 'subdirs' and !$options{Qsub_iterate} );
            next if ( $key eq 'tcga_dirs' and !$options{Qsub_iterate} );
            next if ( !$options{Qsub_iterate} and $options{input_list} and $key eq 'input_list' );    ## If we are in the orignal call, we don't want to qsub more lists
            if ( defined $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}" }
        }
        my $sub_name = $options{sub_name} ? $options{sub_name} : "lgtseqDL";
        ## submit command to grid
        Qsub(
            {   cmd      => $cmd,
                wd       => $lgtseq->{output_dir},
                sub_name => $sub_name,
                sub_mem  => $options{sub_mem},
                sub_mail => $options{sub_mail},
                threads  => $lgtseq->{threads},
                project  => $lgtseq->{project},
            }
        );
        ## Skip to next input for qsub
        my $qsub_pause = defined $options{qsub_pause} ? $options{qsub_pause} : 0;
        sleep $qsub_pause;
        next;
    }
    print_notebook( \%options );

    # Start downloading the input
    if ( defined $lgtseq->{verbose} and $lgtseq->{verbose} == 1 ) { print STDERR "======== lgtseq_download.pl: Starting to download: $input +++\n"; }
    my $downloads;
    ## Process cgquery.xml file
    if ( -e $input and $input =~ /\.xml$/ ) {
        $downloads = $lgtseq->downloadCGHub(
            {   xml        => $input,
                output_dir => $lgtseq->{output_dir},
                threads    => $lgtseq->{threads},
                rate_limit => $lgtseq->{rate_limit},
                cghub_key  => $lgtseq->{cghub_key},
            }
        );
    }
    ## Process TCGA anlysis ID
    elsif ( defined $input and $input =~ /^\w{8}\-\w{4}\-\w{4}\-\w{4}\-\w{12}$/ and $input !~ /\.bam$|f\w{0,3}q$/ ) {
        $downloads = $lgtseq->downloadCGHub(
            {   analysis_id => $input,
                output_dir  => $lgtseq->{output_dir},
                threads     => $lgtseq->{threads},
                rate_limit  => $lgtseq->{rate_limit},
                cghub_key   => $lgtseq->{cghub_key},
            }
        );
    }
    ## Process EGA ID or filename(s)
    elsif ( defined $input and $input =~ /^\d{5,6}(.f\w{0,3}q.gz.gpg)?/ ) {
        $downloads = $lgtseq->downloadEGA(
            {   input      => $input,
                output_dir => $lgtseq->{output_dir},
            }
        );
    }
    else {
        confess "Error: Unable to properly process input.\n";
    }

    # Launch lgtseq_prelim-filter.pl
    if ( ( defined $options{launch_prelim_filter} and $options{launch_prelim_filter} == 1 ) and -e $downloads->[0] ) {
        my $prelim_dir = $options{prelim_dir} ? $options{prelim_dir} : $lgtseq->{output_dir};
        if ( defined $options{prelim_dir} and ( defined $options{tcga_dirs} and $options{tcga_dirs} == 1 ) ) { $prelim_dir = $prelim_dir . "/$tcga_dir\/"; }
        my $prelim_threads  = defined $options{prelim_threads}  ? $options{prelim_threads}  : $lgtseq->{threads};
        my $prelim_sub_mem  = defined $options{prelim_sub_mem}  ? $options{prelim_sub_mem}  : $lgtseq->{sub_mem};
        my $prelim_sort_mem = defined $options{prelim_sort_mem} ? $options{prelim_sort_mem} : $lgtseq->{sort_mem};
        my $lgtseq_prelim_cmd
            = "/home/ksieber/lgtseq/lgtseq_prelim-filter.pl --Qsub=1 --subdirs=0 --output_dir=$prelim_dir --threads=$prelim_threads --sub_mem=$prelim_sub_mem --sort_mem=$prelim_sort_mem";
        ## Bam input
        if ( $downloads->[0] =~ /\.bam/ ) { $lgtseq_prelim_cmd = $lgtseq_prelim_cmd . " --input=$downloads->[0]"; }
        ## Fastq input
        elsif ( $downloads->[0] =~ /.f\w{0,3}q(.gz)?$/ and scalar(@$downloads) == 2 ) { $lgtseq_prelim_cmd = $lgtseq_prelim_cmd . " --input=$downloads->[0]\,$downloads->[1]"; }
        ## Add other potential options for prelim-filter.pl
        my $prelim_name_sort_input = ( defined $options{prelim_name_sort_input} ) ? $options{prelim_name_sort_input} : 0;
        ## If we have a TCGA analysis-id we must prelim_name_sort_input=1
        $prelim_name_sort_input = 1 if ( defined $input and $input =~ /^\w{8}\-\w{4}\-\w{4}\-\w{4}\-\w{12}$/ and $input !~ /\.bam$|f\w{0,3}q$/ );
        if ( defined $prelim_name_sort_input       and $prelim_name_sort_input == 1 )       { $lgtseq_prelim_cmd = $lgtseq_prelim_cmd . " --name_sort_input=1"; }
        if ( defined $options{prelim_delete_input} and $options{prelim_delete_input} == 1 ) { $lgtseq_prelim_cmd = $lgtseq_prelim_cmd . " --delete_input=1"; }
        if ( defined $options{verbose} )  { $lgtseq_prelim_cmd = $lgtseq_prelim_cmd . " --verbose=$options{verbose}"; }     # If the user pass something different than ~/.lgtseek.conf default
        if ( defined $options{sub_mail} ) { $lgtseq_prelim_cmd = $lgtseq_prelim_cmd . " --sub_mail=$options{sub_mail}"; }
        if ( defined $options{prelim_delete_input} and $options{prelim_delete_input} == 1 ) { $lgtseq_prelim_cmd = $lgtseq_prelim_cmd . " --delete_input=1"; }
        ## Add on lgtseq_prelim-filter.pl launching lgtseq_analysis.pl
        if ( defined $options{launch_analysis} and $options{launch_analysis} == 1 ) {
            my $analysis_threads = defined $options{analysis_threads} ? $options{analysis_threads} : $lgtseq->{threads};
            $lgtseq_prelim_cmd = $lgtseq_prelim_cmd . " --launch_analysis=1 --analysis_threads=$analysis_threads";
        }
        if ( defined $lgtseq->{verbose} and $lgtseq->{verbose} == 1 ) { print STDERR "======== lgtseq_download.pl: Starting lgtseq_prelim-filter.pl on: $downloads->[0] +++\n"; }
        Qsub($lgtseq_prelim_cmd);
    }
    elsif ( ( defined $options{launch_analysis} and $options{launch_analysis} == 1 ) and -e $downloads->[0] ) {
        my $analysis_threads  = defined $options{analysis_threads}  ? $options{analysis_threads}  : $lgtseq->{threads};
        my $analysis_sub_mem  = defined $options{analysis_sub_mem}  ? $options{analysis_sub_mem}  : $lgtseq->{sub_mem};
        my $analysis_sort_mem = defined $options{analysis_sort_mem} ? $options{analysis_sort_mem} : $lgtseq->{sort_mem};
        my $lgtseq_analysis_cmd
            = "/home/ksieber/lgtseq/lgtseq_analysis.pl --Qsub=1 --subdirs=0 --prelim_filter=1  --threads=$analysis_threads --sub_mem=$analysis_sub_mem --sort_mem=$analysis_sort_mem";
        if ( $downloads->[0] =~ /\.bam$/ ) {
            $lgtseq_analysis_cmd = $lgtseq_analysis_cmd . " --input=$downloads->[0]";
        }
        elsif ( $downloads->[0] =~ /.f\w{0,3}q(.gz)?$/ and scalar(@$downloads) == 2 ) {
            $lgtseq_analysis_cmd = $lgtseq_analysis_cmd . " --input=$downloads->[0]\,$downloads->[1]";
        }
        if ( defined $options{analysis_dir} ) { $lgtseq_analysis_cmd = $lgtseq_analysis_cmd . " --output_dir=$options{analysis_dir}"; }
        if ( defined $options{analysis_delete_input} and $options{analysis_delete_input} == 1 ) { $lgtseq_analysis_cmd = $lgtseq_analysis_cmd . " --delete_input=1"; }
        Qsub($lgtseq_analysis_cmd);
    }
}

print_complete( \%options );

sub help {
    die "This script will download data, and start further lgtseq analysis. 
         _____________
    ____/Input Options\\__________________________________________________________________________
    Input must be either a:     TCGA analysis_id, cgquery.xml file, EGA id, or EGA filename
    --input=                    Single input to download.
    --input_list=               List with 1 id per line.
         ______________
    ____/Output Options\\_________________________________________________________________________
    --output_dir|o=             Directory for all output. Example: /path/to/{output_dir}/{tcga_dirs}/{subdirs}/ || /path/to/{output_dir}/{subdirs}/
      --tcga_dirs=              <0|1> [0] 1= Maintain TCGA analysis_id directory structure. If the input is an analysis_id, --tcga_dirs=[1].
      --subdirs=                <0|1> [0] 1= Make the sub-dir prefix in output_dir based on input name. If the input is an analysis_id, --subdirs=[1].         
    
    --help|h                    Basic help info.
    --help_full|?               Full  help info.\n\n";
}

sub help_full {
    die "\n$0 : This script will download data.
         _____________
    ____/Input Options\\_____________________________________________________________________________
    Input must be either a:     TCGA analysis_id, cgquery.xml file, EGA id, or EGA filename
    --input=                    Single input to download.
    --input_list=               List with 1 id per line.
         ______________
    ____/Output Options\\____________________________________________________________________________
    --output_dir|o=             Directory for all output. Example: /path/to/{output_dir}/{tcga_dirs}/{subdirs}/ || /path/to/{output_dir}/{subdirs}/
      --tcga_dirs=              <0|1> [0] 1= Maintain TCGA analysis_id directory structure. If the input is an analysis_id, --tcga_dirs=[1].
      --subdirs=                <0|1> [0] 1= Make the sub-dir prefix in output_dir based on input name. If the input is an analysis_id, --subdirs=[1].         
         ___________
    ____/Run Options\\_______________________________________________________________________________  
    --Qsub|q=                   <0|1> [0] 1= Submit each file from --input_list to the grid.
    --Qsub_iterate|Q=           <0|1> [0] 1= Submit job to the grid to iterate over --input_list.
      --threads=                < # > Number of threads to use for download. Only helps TCGA downloads.
      --sub_name=               [lgtseqDL]
      --sub_mail=               <0|1> [0] 1= email \$user\@som.umaryland.edu when job is complete & with stats. Can also specify --sub_mail=specific\@email.foo
    --rate_limit                [~/.lgtseek.config]
    --cghub_key=                [~/.lgtseek.config]
    --qsub_pause=               < # > [0] Number of seconds to pause between each SGE submission. Suggested 30 for EGA ids. 
         _____________
    ____/Prelim Options\\____________________________________________________________________________
    --launch_prelim_filter=     <0|1> [0] 1= Launch lgtseq_prelim-filter.pl after each download is complete. 
      --prelim_name_sort_input= <0|1> [0] 1= Prelim-filter name sort the input (prereq = name sorted). TCGA analysis_id = 1 by default.
      --prelim_dir=             [--output_dir]
      --prelim_threads=         < # > [--threads]
      --prelim_sort_mem=        [~/.lgtseek.config]
      --prelim_sub_mem=         [~/.lgtseek.config]
      --prelim_delete_input=    <0|1> [0] DELETE input file after hg19 aln and/or Prelim-filter
         ________________
    ____/Analysis Options\\___________________________________________________________________________
    --launch_analysis=          <0|1> [0]       1= Start lgtseq_analysis.pl. This can be used with or without lgtseq_prelim-filter.pl.
      --analysis_dir=           [--output_dir]  Only works if {--launch_prelim_filter=0}
      --analysis_threads=       [--threads]     Specify # threads for lgtseq_analysis.pl only.
      --analysis_sort_mem=      [~/.lgtseek.config]
      --analysis_sub_mem=       [~/.lgtseek.config]
      --analysis_delete_input=  <0|1> [0] DELETE input file after hg19 aln and/or Prelim-filter
    __________________________________________________________________________________________________
    --verbose=                  <0|1> [1] 1= Verbose reporting.
    --help|h                    Basic help info.
    --help_full|?               Full  help info.
    __________________________________________________________________________________________________\n";
}
