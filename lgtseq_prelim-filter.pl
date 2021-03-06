#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/

=head1 NAME

lgtseq_prelim-filter.pl

=head1 SYNOPSIS

Filter a bam before running lgtseq_analysis.pl 

=head1 DESCRIPTION

Prelimary filter out Human reads and split bams into smaller chunks of potetional bacterial DNA integration and microbiome reads before running lgtseq_analysis.pl.
This allow for lgtseq_analysis.pl to run dramatically quicker.

=head1 EXAMPLE

lgtseq_prelim-filter.pl --input=/path/to/unfiltered.bam --output_dir=/dir/for/output/ --qsub=1 --threads=4

=head1 AUTHOR - Karsten Sieber

e-mail: Karsten.Sieber@gmail.com

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

my $LGTSEQ_PRELIM = '1.07';

use warnings;
no warnings 'uninitialized';
use strict;
use Carp;
$Carp::MaxArgLen = 0;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
our %options;
my $results = GetOptions(
    \%options,         'input_list|I=s',  'input|i=s',         'name_sort_input=s',  'sort_mem=s',          'split_bam=i',   'encrypt=i',        'key=s',
    'prelim_filter=i', 'seqs_per_file=i', 'keep_softclip=i',   'Qsub|q=i',           'excl=i',              'sub_mem=s',     'sub_name=s',       'threads|t=i',
    'projects=s',      'output_dir|o=s',  'subdirs=i',         'overwrite=s',        'samtools_bin=s',      'ergatis_dir=s', 'output_list=s',    'bin_dir=s',
    'fs=s',            'clovr=s',         'diag',              'verbose|V=i',        'print_hostname|ph=i', 'config_file=s', 'help|h',           'help_full|?',
    'tcga_dirs=i',     'aln_human=i',     'hg19_ref=s',        'no_gal=i',           'hostname=s',          'sub_mail=s',    'Qsub_iterate|Q=i', 'launch_analysis=i',
    'analysis_dir=s',  'analysis_iter=i', 'name_sort_check=i', 'analysis_threads=i', 'analysis_mem=s',      'qsub_pause=i',  'delete_input=i',   'analysis_subdirs=i',
) or die "Error: Unrecognized command line option. Please try again.\n";
use print_call;

# print_hostname(\%options);                                    ## This is useful for trouble shooting grid nodes that might be missing modules etc.
use Scalar::Util qw(reftype);
use POSIX;
use run_cmd;
use lib ( '/home/ksieber/perl5/lib/perl5/', '/local/projects-t3/HLGT/scripts/lgtseek/lib/', '/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/' );
use LGTSeek;
use File::Basename;
use setup_input;
use setup_output;

if ( $options{help} )      { &help; }
if ( $options{help_full} ) { &help_full; }

if ( !$options{input} and !$options{input_list} ) { die "Error: Must give input with --input=<BAM> or --input_list=<LIST of bams>\n"; }

## Set default values
$options{prelim_filter} = defined $options{prelim_filter} ? "$options{prelim_filter}" : "1";
$options{sub_mem}       = defined $options{sub_mem}       ? "$options{sub_mem}"       : "5G";
$options{output_list}   = defined $options{output_list}   ? "$options{output_list}"   : "1";

if ( defined $options{input_list} ) { $options{subdirs} = 1; }

my $lgtseek = LGTSeek->new2( \%options );
my $delete_input = defined $options{delete_input} ? $options{delete_input} : $lgtseek->{delete_input};  ## Currently set to 0

my $inputs;
if ( !$options{Qsub_iterate} ) {
    $inputs = setup_input( \%options );
}
elsif ( $options{input_list} ) {
    push( my @tmp_input_list, $options{input_list} );
    $inputs = \@tmp_input_list;
}

my $x = 0;

my $original_output_dir = $lgtseek->{output_dir};
my $fastqs;

## Print the script call
print_call( \%options, "LGTSEQ_PRELIM_VERSION=$LGTSEQ_PRELIM\tLGTSeek.pm_VERSION=$LGTSeek::VERSION" );

foreach my $input (@$inputs) {

    # Set a few defaults
    ## Defaults for TCGA analysis-id input
    if ( $input =~ /^\w{8}\-\w{4}\-\w{4}\-\w{4}\-\w{12}$/ and ( $input !~ /\.bam$/ and $input !~ /\.f\w{0,3}q(.gz)?$/ ) ) {
        $lgtseek->{tcga_dirs} = 1;        # If we have a TCGA analysis ID, make sure --tcga_dirs is on.
        $options{name_sort_input} = 1;    # We know if we download the file we have to name sort it.
    }
    ## Defaults for EGA-id input and split comma seperated fastqs
    elsif ( ( $input =~ /^\d{5,6}.*(\.f\w{0,3}q\.gz)?$/ or $input =~ /(\.f\w{0,3}q\.gz)?$/ ) and $input !~ /\.bam$/ ) {
        my ( $in1, $in2 ) = split( /,/, $input );
        $fastqs             = $input;     ## Keep potentially comma seperated _1 & _2 for later
        $input              = $in1;       ## Change input to _1 for clean File::Basename parsing
        $options{aln_human} = 1;
    }

    my ( $fn, $path, $suf ) = fileparse( $input, ( @{ $lgtseek->{bam_suffix_list} }, @{ $lgtseek->{fastq_suffix_list} } ) );
    my $subdir = ( $suf =~ /\w+/ ) ? $fn : undef;
    my @split_path = split( /\//, $path );
    my $tcga_dir = ( $suf =~ /\.bam$|\.f\w{0,3}q(.gz)?/ ) ? $split_path[-1] : $fn;
    $lgtseek->{output_dir} = defined $original_output_dir ? $original_output_dir : $path;

    if ( $lgtseek->{tcga_dirs} == 1 and !$options{Qsub_iterate} ) {
        $lgtseek->{output_dir} = $lgtseek->{output_dir} . "/$tcga_dir\/" unless ( $lgtseek->{output_dir} =~ /$tcga_dir\/*$/ );
    }

    if ( $lgtseek->{subdirs} == 1 and !$options{Qsub_iterate} and defined $subdir ) {
        $lgtseek->{output_dir} = $lgtseek->{output_dir} . "/$subdir\/" unless ( $lgtseek->{output_dir} =~ /$subdir\/*$/ );
    }
    $lgtseek->{output_dir} =~ s/\/{2,}/\//g;
    $options{output_dir} = $lgtseek->{output_dir};
    $lgtseek->_run_cmd("mkdir -p -m u=rwx,g=rwx,o= $lgtseek->{output_dir}");
    chdir( $lgtseek->{output_dir} );

    ## Qsub this script foreach input and any of the options passed
    if ( $lgtseek->{Qsub} == 1 or $options{Qsub_iterate} == 1 ) {
        ## If we are in the orignal call, change input from list to a single file
        if ( $options{input_list} and !$options{Qsub_iterate} ) { $options{input} = $input; }
        ## Check $sub_mem is enough for sorting
        if ( $lgtseek->{name_sort_input} == 1 ) {
            my $original_sub_mem;
            my $original_sort_mem;
            if ( $lgtseek->{sub_mem} =~ /^(\d+)[K|M|G]$/ )  { $original_sub_mem  = $1; }
            if ( $lgtseek->{sort_mem} =~ /^(\d+)[K|M|G]$/ ) { $original_sort_mem = $1; }
            if ( $original_sub_mem < ( $original_sort_mem * $options{threads} ) ) {
                $options{sub_mem} = ( ceil( ( $original_sort_mem * $options{threads} ) * 1.1 ) ) + 1 . "G";
            }
        }

        ## Build qsub command
        my $cmd = "$^X $0";
        foreach my $key ( keys %options ) {
            next if ( $key eq 'Qsub' );
            next if ( $key eq 'subdirs' and !$options{Qsub_iterate} );
            next if ( $key eq 'tcga_dirs' and !$options{Qsub_iterate} );
            next if ( $key eq 'Qsub_iterate' );
            next if ( !$options{Qsub_iterate} and $options{input_list} and $key eq 'input_list' );    ## If we are in the orignal call, we don't want to qsub more lists
            if ( defined $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}" }
        }
        my $sub_name   = $options{sub_name}           ? $options{sub_name}   : "prelim$x";
        my $qsub_pause = defined $options{qsub_pause} ? $options{qsub_pause} : 0;
        ## submit command to grid
        Qsub(
            {   cmd      => $cmd,
                wd       => $lgtseek->{output_dir},
                sub_name => $sub_name,
                sub_mem  => $lgtseek->{sub_mem},
                sub_mail => $options{sub_mail},
                threads  => $lgtseek->{threads},
                project  => $lgtseek->{project},
                hostname => $options{hostname},
            }
        );
        sleep $qsub_pause;
        $x++;
        ## Skip to next input for qsub
        next;
    }
    if ( scalar( @{$inputs} ) > 1 ) { print_call( \%options ); }
    print_notebook( \%options );

    my $bam_downloaded;
    if ( $input =~ /^\w{8}\-\w{4}\-\w{4}\-\w{4}\-\w{12}$/ and $input !~ /\.bam$|f\w{q$/ ) {
        $bam_downloaded = $lgtseek->downloadCGHub(
            {   analysis_id => $input,
                output_dir  => $lgtseek->{output_dir},
                threads     => $lgtseek->{threads},
                rate_limit  => $lgtseek->{rate_limit},
                cghub_key   => $lgtseek->{cghub_key},
            }
        );
        $input = $bam_downloaded->[0];
    }
    elsif ( !-e $input and ( $input =~ /^\d{5,6}(.f\w{0,3}q.gz.gpg)?$/ and $input !~ /\.bam$/ ) ) {
        $bam_downloaded = $lgtseek->downloadEGA(
            {   input      => $input,
                output_dir => $lgtseek->{output_dir},
            }
        );
        $input = $bam_downloaded->[0];
    }

    # Align to Human.
    my $human_aln = defined $options{aln_human} ? $options{aln_human} : "0";
    if ( $human_aln == 1 ) {
        print STDERR "======== RUNBWA-Human ========\n";
        my $human_bam1;
        ## Map input bam @ hg19
        if ( $input =~ /\.bam$/ ) {
            $human_bam1 = $lgtseek->runBWA(
                {    ## &runBWA returns an array
                    input_bam   => $input,
                    output_bam  => 1,
                    threads     => $lgtseek->{threads},
                    output_dir  => "$lgtseek->{output_dir}/aln_human/",
                    reference   => $lgtseek->{hg19_ref},
                    overwrite   => $lgtseek->{overwrite},
                    cleanup_sai => 1,
                }
            );
            if ( defined $delete_input and $delete_input == 1 ) {
                $lgtseek->_run_cmd("rm $input") if ( -e $input );
            }
        }
        ## Map fastqs @ hg19
        elsif ( $input =~ /\.f\w{0,3}q(.gz)?$/ ) {
            my ( $in1, $in2 ) = split( /,/, $fastqs );    ## split $fastqs input if needed
            my ( $input_base, $input_dir, $input_suffix ) = fileparse( $in1, @{ $lgtseek->{fastq_suffix_list} } );
            ## &runBWA returns an array
            $human_bam1 = $lgtseek->runBWA(
                {   input_dir   => $input_dir,
                    input_base  => $input_base,
                    output_bam  => 1,
                    threads     => $lgtseek->{threads},
                    output_dir  => "$lgtseek->{output_dir}/aln_human/",
                    reference   => $lgtseek->{hg19_ref},
                    overwrite   => $lgtseek->{overwrite},
                    cleanup_sai => 1,
                }
            );
            if ( defined $delete_input and $delete_input == 1 ) {
                $lgtseek->_run_cmd("rm $in1") if ( -e "$in1" );
                $lgtseek->_run_cmd("rm $in2") if ( -e "$in2" );
            }
        }
        $input = @$human_bam1[0];
    }

    ## PrelimFilter to remove M_M reads.
    my $bams;
    if ( $lgtseek->{prelim_filter} == 1 || $lgtseek->{split_bam} == 1 || $lgtseek->{name_sort_input} == 1 ) {
        $bams = $lgtseek->prelim_filter(
            {   input_bam       => $input,
                output_dir      => $lgtseek->{output_dir},
                name_sort_input => $lgtseek->{name_sort_input},    ## Default = 0
                name_sort_check => $options{name_sort_check},      ## Default = $name_sort_input
                sort_mem        => $lgtseek->{sort_mem},           ## Default = 1G lgtseek default.
                split_bam       => $lgtseek->{split_bam},          ## Default = 1
                seqs_per_file   => $lgtseek->{seqs_per_file},      ## Default = 50M
                keep_softclip   => $lgtseek->{keep_softclip},      ## Default = 1
                overwrite       => $lgtseek->{overwrite},          ## Default = 0
            }
        );
        if ( defined $delete_input and $delete_input == 1 ) {
            my ( $input_fn, $input_path, $input_suffix ) = fileparse( $input, ".bam" );
            $lgtseek->_run_cmd("rm $input")      if ( -e "$input" );
            $lgtseek->_run_cmd("rm $input\.bai") if ( -e "$input\.bai" );
            if ( -e $input_path and $input_path =~ /$tcga_dir\/$tcga_dir/ ) { $lgtseek->_run_cmd("rmdir --ignore-fail-on-non-empty $input_path"); }
            if ( $human_aln == 1 ) { $lgtseek->_run_cmd("rm -rf $lgtseek->{output_dir}/aln_human/"); }
        }

    }

    ## Encrypt the file(s)
    my @encrypted;
    if ( $lgtseek->{encrypt} == 1 ) {
        foreach my $files (@$bams) {
            my $cmd = "gpg -o $files\.gpg --default-key $options{key} -r $options{key} -e $files";
            if   ( $lgtseek->{Qsub} == 1 ) { Qsub($cmd); }
            else                           { $lgtseek->_run_cmd($cmd); }
            push( @encrypted, "$files\.gpg" );
        }
    }

    ## Print out a list of the _prelim.bams that were created
    if ( $lgtseek->{output_list} == 1 ) {
        open( my $olistfh, ">$lgtseek->{output_dir}/output.list" ) or confess "Unable to open: $lgtseek->{output_dir}/output.list because: $!\n";
        if ( $lgtseek->{encrypt} == 1 ) {
            foreach my $out (@encrypted) { print $olistfh "$out\n"; }
        }
        elsif ( $lgtseek->{split_bam} == 1 || $lgtseek->{prelim_filter} == 1 ) {
            foreach my $out2 (@$bams) { print $olistfh "$out2\n"; }
        }
    }

    if ( defined $options{launch_analysis} and $options{launch_analysis} == 1 ) {
        my $analysis_dir = $options{analysis_dir} ? $options{analysis_dir} : $lgtseek->{output_dir};
        if ( defined $options{analysis_dir} and ( defined $options{tcga_dirs} and $options{tcga_dirs} == 1 ) ) { $analysis_dir = $analysis_dir . "/$tcga_dir\/"; }
        my $analysis_threads = defined $options{analysis_threads} ? $options{analysis_threads} : $lgtseek->{threads};
        my $analysis_mem     = defined $options{analysis_mem}     ? $options{analysis_mem}     : "5G";
        my $analysis_subdirs = defined $options{analysis_subdirs} ? $options{analysis_subdirs} : "1";
        my $lgtseq_analysis_cmd
            = "/home/ksieber/scripts/lgtseq_analysis.pl --input_list=$lgtseek->{output_dir}/output.list --output_dir=$analysis_dir --threads=$analysis_threads --sub_mem=$analysis_mem --subdirs=$analysis_subdirs";
        if   ( $options{analysis_iter} == 1 ) { $lgtseq_analysis_cmd = $lgtseq_analysis_cmd . " --Qsub_iterate=1"; }
        else                                  { $lgtseq_analysis_cmd = $lgtseq_analysis_cmd . " --Qsub=1"; }
        if ( $options{verbose} == 1 )  { $lgtseq_analysis_cmd = $lgtseq_analysis_cmd . " --verbose=1"; }
        if ( $options{sub_mail} == 1 ) { $lgtseq_analysis_cmd = $lgtseq_analysis_cmd . " --sub_mail=$options{sub_mail}"; }
        if ( $options{verbose} ) { print STDERR "======== &prelim_filter: Starting lgtseq_analysis.pl on: $lgtseek->{output_dir}/output.list. +++\n"; }
        Qsub( { cmd => $lgtseq_analysis_cmd } );
    }
}

&print_complete( \%options, "LGTSEQ_PRELIM_VERSION=$LGTSEQ_PRELIM\tLGTSeek.pm_VERSION=$LGTSeek::VERSION" );

sub help {
    die "Help Basic Info: This script will remove M_M reads and split output bams into smaller bams. 
        --input=            <BAM>
        --name_sort_input=  <0|1> [0] 1= Resort the input bam by read names.  
        --output_dir=       Directory for output. 
        --help_full|?       Full Help Info\n";
}

sub help_full {
    die "   LGTSEQ_PRELIM_VERSION  $LGTSEQ_PRELIM
    Help Full Info: This script will remove M_M reads, keeping M_U, U_U, M_M with Softclip. The output bams are split into smaller chunks.
    This script is primarily used to filter data before lgtseq_analysis.pl. It can also download CGHub data, filter, and start lgtseq_analysis.pl
             _____________
        ____/Input Options\\__________________________________________________________________________
        --input|i=              File to LGTSeq prelim filter. May be: .bam, .fq, or an TCGA analysis_id to download. 
        --input_list|I=         <LIST of INPUTS> File with 1 input per line or comma seperated at command line.
        --delete_input=         <0|1> [1] 1= DELETE input after alignment to hg19 and/or prelim_filtering.
             _________________
        ____/Filtering Options\\______________________________________________________________________
        --aln_human =           <0|1> [0] 1= Align the input bam or fastq's to hg19 before filtering. **MUST** be name sorted input.
        --prelim_filter=        <0|1> [1] 1= Filter out M_M reads, keeping MU,UU,and SC. 
          --keep_softclip=      <0|1> [1] 1= Keep soft clipped reads >=24 bp (Reads potentially on LGT) 
        --name_sort_input=      <0|1> [0] 1= Resort the input bam by read names. If the input is analysis_id, --name_sort_input=[1].
          --sort_mem=           [1G] Mem per thread to sort with. 
          --threads=            [1] # of threads. 
        --name_sort_check=      <0|1> [--name_sort_input] 1= Quick & dirty check for proper pairing of sorted.bam. 
        --split_bam=            <0|1> [1] 1= Split bam(s)
          --seqs_per_file=      <lines per bam> [50,000,000]
             ______________
        ____/Output Options\\_________________________________________________________________________
        --output_dir|o=         Directory for all output. Example: /path/to/{output_dir}/{tcga_dirs}/{subdirs}/ || /path/to/{output_dir}/{subdirs}/
         --tcga_dirs=           <0|1> [0] 1= Maintain TCGA analysis_id directory structure. If the input is an analysis_id, --tcga_dirs=[1].
          --subdirs=            <0|1> [0] 1= Make the sub-dir prefix in output_dir based on input name. If the input is an analysis_id, & --subdirs=[1].
        --overwrite=            <0|1> [1] 1= Overwrite previous files. ** HIGHLY ** Recommended to =1 (ON) if *restarting* lgtseq_prelim-filter.pl.
             ___________
        ____/SGE Options\\____________________________________________________________________________
        --Qsub|q=               <0|1> [0] 1= Qsub this script for each input.
        --Qsub_iterate=         <0|1> [0] 1= Submit job to the grid to iterate over --input_list.
          --project=            <project> [jdhotopp-lab] SGE group project to submit command under.
          --sub_mem=            [5G] --sub_mem MUST >= (--threads * --sort_mem)
          --sub_name=           < > Name of the SGE Job.
          --sub_mail=           [0] 1= email \$user\@som.umaryland.edu when job is complete & with stats. Can also specify --sub_mail=specific\@email.foo
          --qsub_pause=         < # > [0] Number of seconds to pause between each SGE submission. 
             _____________
        ____/Launch LGTSeq\\__________________________________________________________________________
        --launch_analysis=      <0|1> [0]       1= Start lgtseq_analysis.pl on each file from the prelim_filtering output.list.  
                                  ~/.lgtseek.conf defaults, except {--threads}, {--sub_mem}, {--subdirs=1}, & {--analysis_dir}
        --analysis_dir=         [--output_dir]  Specify the directory for lgtseq_analysis.
        --analysis_threads=     [--threads]     Specify # threads for lgtseq_analysis.pl only.
        --analysis_mem=         [.lgtseq.conf]  Specify SGE sub mem for lgtseq_analysis.pl. 
        --analysis_iter=        <0|1> [0]       1= The lgtseq_analysis.pl will iterate over the list of outputs instead of submitting each seperately. 
                                                0= Submit each file from the output.list from prelim_filtering to the grid for lgtseq_analysis. 
                                                Suggested NOT to --split_bams if iterating. 
        ----------------------------------------------------------------------------------------------
        --verbose               <0|1> [1] 1= Verbose reporting of progress. 0 =Turns off reports. 
        --help|h                Basic Help Info
        --help_full|?           Full Help Info
        --conf_file=            [~/.lgtseek.conf]
        ----------------------------------------------------------------------------------------------\n";
}
__END__
