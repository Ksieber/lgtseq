#!/usr/bin/perl -I /home/ksieber/perl5/lib/perl5/ -I /home/ksieber/scripts/

my $VERSION = "1.01";

use warnings;
no warnings 'uninitialized';
use strict;
use lib ( "/local/projects-t3/HLGT/scripts/lgtseek/lib/", "/local/projects/ergatis/package-driley/lib/perl5/x86_64-linux-thread-multi/" );
use LGTSeek;
use run_cmd;
use mk_dir;
use print_call;
use setup_input;
use Time::SoFar;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
my %options;
my $results = GetOptions(
    \%options, 'input|i=s',  'input_list|I=s', 'output_dir|o=s', 'output_prefix|p=s', 'Qsub|q=i',      'subdirs=s',              'conf_file=s',
    'help|?',  'sub_name=s', 'sub_mail=s',     'iter=i',         'merge=i',           'append_fn|F=i', 'append_analysis_id|A=i', 'append_tcga_id|T=i',
);

if ( $options{help} ) {
    die "Help:
    This script will merge the output files of lgt_seq.pl.
    ----------------------------------------------------------------------------------------
    --input|i=                  Directory to merge files from. Subdirectories in this directory should contain lgtseq_analysis.pl outputs.
    --input_list|I=             List of input directories.
     --merge=                   <0|1> [1] 1 = Will merge all of the directories in the list into one.
     --iter=                    <0|1> [0] 1 = Will iterate over all the directories in the list. 
    ----------------------------------------------------------------------------------------
    --output_dir|o=             Directory for output.
      --subdirs=                <0|1> [0] 1= Make a new directory within output_dir for each input directory.
      --output_prefix|p=        Prefix name for the merged output. [Input dir name]
      --append_fn|F=            <0|1> [1] Append file each read originated from. 
      --append_analysis_id|A=   <0|1> [1] Append analysis-id to each read.    \$path =~ m/(\\w{8}-\\w{4}-\\w{4}-\\w{4}-\\w{12})/;
      --append_tcga_id|T=       <0|1> [1] Append tcga-id to each read.        \$fn   =~ m/(TCGA-\\w{2}-\\w{4}-\\w{3}-\\w{3}-\\w{4})/;
    ----------------------------------------------------------------------------------------
    --Qsub|q=i                  <0|1> [0] 1= Qsub the job to the grid. 
     --sub_mail=                [0] 1= email user\@som.umaryland.edu when job is complete & with stats. Can also specify --sub_mail=specific\@email.foo
    ----------------------------------------------------------------------------------------
    --help|?
    --conf_file=                [~/.lgtseek.conf]
    ----------------------------------------------------------------------------------------\n";
}

if ( !$options{input} && !$options{input_list} ) { die "Must give an input. Use --input or --input_list\n"; }
if ( !$options{output_dir} ) { die "Must use --output_dir=\n"; }
my $script_command_line = build_call_string( \%options );
my $iter = defined $options{iter} ? $options{iter} : "0";
my $merge;
if    ( $iter && !$options{merge} ) { $merge = "0"; }
elsif ( $iter && $options{merge} )  { $merge = $options{merge}; }
else                                { $merge = "1"; }

my $subdirs = defined $options{subdirs} ? $options{subdirs} : "0";
if ( $options{input_list} && $iter ) { $subdirs = 1; }
$options{append_fn}          = defined $options{append_fn}          ? $options{append_fn}          : 1;
$options{append_analysis_id} = defined $options{append_analysis_id} ? $options{append_analysis_id} : 1;
$options{append_tcga_id}     = defined $options{append_tcga_id}     ? $options{append_tcga_id}     : 1;

my $input = setup_input( \%options );
run_cmd("mkdir -p $options{output_dir}");

if ( $options{Qsub} ) {
    my $sub_name = defined $options{sub_name} ? $options{sub_name} : "mergelgtseq";
    my $cmd = "$^X $0";
    foreach my $key ( keys %options ) {
        next if ( $key eq 'Qsub' );
        next if ( $key eq 'sub_name' );
        if ( $options{$key} ) { $cmd = $cmd . " --$key=$options{$key}" }
    }
    Qsub( { cmd => $cmd, sub_name => $sub_name, sub_mail => $options{sub_mail} } );
    die "+++ Job submitted to the grid. +++\n";
}

my $lgtseek = LGTSeek->new2( \%options );

my @bam_suffix_to_merge = (
    'microbiome.bam',                  'lgt_host.bam',                     'integration_site_donor_host.bam', 'lgt_host_filtered.bam',
    'lgt_host_filtered_validated.bam', 'integration_site_donor_donor.bam', 'microbiome_filtered.bam',         'lgt_final.bam'
);

my @txt_suffix_to_merge = (
    'by_clone.txt',          'by_trace.txt',                           'post_processing.tab',    'lgt_host_lineage1.out',
    'lgt_host_lineage2.out', 'microbiome_lca-bwa_independent_lca.txt', 'microbiome_lca-bwa.txt', 'lgt_lca-bwa.txt',
    'lgt__lca-bwa.txt'
);

if ($merge) { &merge_dirs($input); }
if ($iter)  { &iter_dirs($input); }

print_complete( \%options );

###################
### Subroutines ###
###################

sub iter_dirs {
    my $input = shift;
    foreach my $dir (@$input) {
        my @split_dir  = split( /\//, $dir );
        my $name       = $split_dir[-1];
        my $output_dir = "$options{output_dir}\/$name";
        if ( -e $output_dir ) { $output_dir = "$options{output_dir}\/$name/merged/"; }
        mk_dir($output_dir);

        ## Process all the txt files for merging
        # print STDERR "Input_dir: $dir\noutput_prefix: $name\nOutput_dir: $output_dir\n";
        foreach my $txt_suffix (@txt_suffix_to_merge) {
            chomp( my @list_to_merge = `find $dir -name '*$txt_suffix'` );
            if ( scalar( @list_to_merge == 0 ) ) { print STDERR "*** Warning *** : Did not find any *$txt_suffix in: $dir\n"; next; }
            my $output = "$output_dir/$name\_$txt_suffix";
            if ( $txt_suffix =~ /by_trace.txt/ ) {
                run_cmd("head -n1 $list_to_merge[0] > $output");    # Grab the header and put it into the output file
                foreach my $file (@list_to_merge) {
                    next if ( run_cmd("head -n2 $file | wc -l") != 2 );
                    run_cmd("grep -v Read $file >> $output");       # Merge all the files together skipping the header
                }
            }
            else {
                foreach my $file (@list_to_merge) {
                    run_cmd("cat $file >> $output");
                }
            }
        }

        ## Process all the bam files for merging
        foreach my $bam_suffix (@bam_suffix_to_merge) {
            chomp( my @list_to_merge = `find $dir -name '*$bam_suffix'` );
            if ( scalar( @list_to_merge == 0 ) ) { print STDERR "*** Warning *** : Did not find any *$bam_suffix in: $dir\n"; next; }
            my $output = "$output_dir/$name\_$bam_suffix";
            my $header;    ## Trying to come up with a way to grab a header while checking to make sure at least 1 file has data in it.
            for ( my $i = 0; $i < scalar @list_to_merge; $i++ ) {
                next if ( $lgtseek->empty_chk( { input => $list_to_merge[$i] } ) == 1 );    ## Skip the bam if it is empty
                $header = merge_headers( $header, $list_to_merge[$i] );
            }
            next if ( !$header->{'SQ'} );                                                   ## If none of the bams in the @list_to_merge have data header should be undef still and we skip this suffix
            open( my $out, "| samtools view -S - -bo $output" ) or die "Can not open output: $output\n";
            print_mereged_header( $header, $out );
            foreach my $bam (@list_to_merge) { process_bam( $bam, $out ); }
            close $out or die "Can't close output: $output because: $!\n";
        }
        print STDERR "====== Completed merging: $dir ======\n";
    }
}

sub merge_dirs {
    my $input_path;
    if ( defined $options{input_list} ) {
        my ( $fn, $path, $suf ) = fileparse( $options{input_list} );
        my @split_path = split( /\//, $path );
        $input_path = $split_path[-1];
    }
    elsif ( $options{input} && -e $options{input} ) {
        my @split_path = split( /\//, $options{input} );
        $input_path = $split_path[-1];
    }
    my $name = defined $options{output_prefix} ? $options{output_prefix} : $input_path;
    my $output_dir;
    if ( $subdirs == 1 ) {
        $output_dir = "$options{output_dir}/$name";
        if ( -e $output_dir ) { $output_dir = "$options{output_dir}/$name/merged/"; }
    }
    else {
        $output_dir = "$options{output_dir}";
    }

    mk_dir($output_dir);

    my $input = shift;

    foreach my $txt_suffix (@txt_suffix_to_merge) {
        my @list_to_merge;
        foreach my $dir (@$input) {
            chomp( my @tmp_list_to_merge = `find $dir -name '*$txt_suffix'` );
            if ( scalar( @tmp_list_to_merge == 0 ) ) { print STDERR "* Warning * : Did not find any \*$txt_suffix in :$dir\n"; next; }
            push( @list_to_merge, @tmp_list_to_merge );
        }
        if ( scalar( @list_to_merge == 0 ) ) { print STDERR "*** Warning *** : Did not find any \*$txt_suffix\n"; next; }
        my $output = "$output_dir/$name\_$txt_suffix";
        if ( $txt_suffix =~ /by_trace.txt/ ) {
            run_cmd("head -n1 $list_to_merge[0] > $output");    # Grab the header and put it into the output file
            foreach my $file (@list_to_merge) {
                next if ( run_cmd("head -n2 $file | wc -l") != 2 );
                run_cmd("grep -v Read $file >> $output");       # Merge all the files together skipping the header
            }
        }
        else {
            foreach my $file (@list_to_merge) {
                run_cmd("cat $file >> $output");
            }
        }
    }

    ## Process all the bam files for merging
    foreach my $bam_suffix (@bam_suffix_to_merge) {
        my %bams_to_merge;
        foreach my $dir (@$input) {
            chomp( my @tmp_list_to_merge = `find $dir -name '*$bam_suffix'` );
            if ( scalar( @tmp_list_to_merge == 0 ) ) { print STDERR "* Warning * : Did not find any *$bam_suffix in: $dir\n"; next; }
            foreach my $file (@tmp_list_to_merge) { $bams_to_merge{$file}++; }
        }
        if ( !%bams_to_merge ) { print STDERR "*** Warning *** : Did not find any *$bam_suffix\n"; next; }
        my $output = "$output_dir/$name\_$bam_suffix";
        my $header;    ## Trying to come up with a way to grab a header while checking to make sure at least 1 file has data in it.
        foreach my $bam_to_check ( keys %bams_to_merge ) {
            next if ( $lgtseek->empty_chk( { input => $bam_to_check } ) == 1 );    ## Skip the bam if it is empty
            $header = &merge_headers( $header, $bam_to_check );
        }
        next if ( !$header->{'SQ'} );                                              ## If none of the bams in the @list_to_merge have data header should be undef still and we skip this suffix
        open( my $out, "| samtools view -S - -bo $output" ) or die "Can not open output: $output\n";
        &print_mereged_header( $header, $out );
        foreach my $bam ( keys %bams_to_merge ) { &process_bam( $bam, $out ); }
        close $out or die "Can't close output: $output because: $!\n";
    }
    print STDERR "====== Completed LGTSEQ-Merging ======\n";

}

sub merge_headers {
    my $header_obj = shift;
    my $file       = shift;
    open( my $header_fh, "-|", "samtools view -H $file" ) or die "Error: Unable to the header on this bam: $file\n";
    my $counter = 0;
    while (<$header_fh>) {
        chomp( my $line = $_ );
        my @split_line = split( /\s+/, $line );
        next if ( $split_line[0] =~ /\@HD/ );
        next if ( $split_line[0] =~ /\@PG/ );
        if ( $split_line[0] =~ /\@SQ/ ) {
            $split_line[1] =~ /SN:([A-Za-z0-9\.\-\_\|]+)/;
            my $chr = $1;
            $header_obj->{'SQ'}->{$chr} = $line;
        }
        if ( $split_line[0] =~ /\@CO/ ) {
            $line =~ s/\s+(PP:[A-Za-z0-9\-\_\.]+[ A-Za-z0-9\-\_\.]?)//;
            if ( !$header_obj->{'CO'}->{$line} ) { $counter++; $header_obj->{'CO'}->{$counter} = $line }
        }
    }
    close $header_fh;
    return $header_obj;
}

sub print_mereged_header {
    my $header_obj = shift;
    my $output_fh  = shift;

    print $output_fh "\@HD\tVN:$LGTSeek::VERSION\tSO:unsorted\n";
    foreach my $chr (
        sort { ( ( $a =~ /[chr]?(\d+)\_?[A-Za-z\-]*(\d*)[A-Za-z\_\-]*/ )[0] or 999 ) <=> ( ( $b =~ /[chr]?(\d+)\_?[A-Za-z\-]*(\d*)[A-Za-z\_\-]*/ )[0] or 999 ) }
        keys %{ $header_obj->{'SQ'} }
        )
    {
        print $output_fh "$header_obj->{'SQ'}->{$chr}\n";
    }
    foreach my $co ( sort keys %{ $header_obj->{'CO'} } ) { print $output_fh "$header_obj->{'CO'}->{$co}\n"; }
    print $output_fh "\@CO\tID:Merged\tPG:lgtseq_merge.pl\tVN:$VERSION\tCL:$script_command_line\n";

}

sub process_bam {
    my $bam    = shift;
    my $out_fh = shift;

    return if ( $lgtseek->empty_chk( { input => $bam } ) == 1 );
    open( my $in_fh, "-|", "samtools view $bam" ) or die "Can not open input: $bam\n";
    my ( $bam_fn, $bam_path, $bam_suf ) = fileparse( $bam, @{ $lgtseek->{bam_suffix_list} } );
    while (<$in_fh>) {
        chomp( my $bam_line = $_ );
        print $out_fh "$bam_line";
        if ( defined $options{append_fn}          and $options{append_fn} == 1 )          { print $out_fh "\tFN:Z:$bam"; }
        if ( defined $options{append_tcga_id}     and $options{append_tcga_id} == 1 )     { print $out_fh "\tTI:Z:$1" if ( $bam_fn =~ /(TCGA\-\w{2}\-\w{4}\-\w{3}\-\w{3}\-\w{4})/ ); }
        if ( defined $options{append_analysis_id} and $options{append_analysis_id} == 1 ) { print $out_fh "\tAI:Z:$1" if ( $bam_path =~ /(\w{8}\-\w{4}\-\w{4}\-\w{4}\-\w{12})/ ); }
        print $out_fh "\n";
    }
    close $in_fh or die "Can not close input: $input because: $!\n";
}

sub build_call_string {
    my $options = shift;
    my $call    = "$0";
    if   ( defined $options->{input} ) { $call = $call . " --input=$options->{input}"; }
    else                               { $call = $call . " --input_list=$options->{input_list}"; }
    foreach my $key (qw( merge iter output_dir subdirs output_prefix append_fn append_analysis_id append_tcga_id Qsub sub_mail)) {
        if ( defined $options->{$key} ) { $call = $call . " --$key=$options->{$key}"; }
    }
    return $call;
}

