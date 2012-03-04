#!/usr/bin/env perl

use strict;
use warnings;
use lib map{'/opt/local/lib/perl5/site_perl/5.8.9'.$_} ('','/darwin-2level');
use DBI;
use CGI qw(:standard escape unescape);
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use List::Util qw(sum);

my $dbh_arabidopsis = Gramene::DB->new('diversity_arabidopsis');
my $dbh_ontology 	= Gramene::DB->new('ontology');

my $query_arabidopsis = <<QUERY;
    SELECT
        local_trait_name,
        to_accession
    FROM
        div_trait_uom
    GROUP BY
        local_trait_name
QUERY

my $query_ontology = <<QUERY;
    SELECT
        term_name
    FROM
        term
    WHERE
        term_accession=?
QUERY

my $arrayref_arabidopsis = $dbh_arabidopsis->selectall_arrayref($query_arabidopsis, { Columns => {} });
my @array_arabidopsis    = @$arrayref_arabidopsis;
my %traits;

for ( my $i = 0; $i < @array_arabidopsis; $i++ ) {

	my $to       = $array_arabidopsis[$i]->{'to_accession'};
	my $trait    = $array_arabidopsis[$i]->{'local_trait_name'};

	$to =~ s/^\s+//;
	$to =~ s/\s+$//;
	$trait =~ s/^\s+//;
	$trait =~ s/\s+$//;
	
	my $arrayref_ontology = $dbh_ontology->selectall_arrayref($query_ontology, { Columns => {} }, ( $to ));
	my @array_ontology    = @$arrayref_ontology;
	
	my $name  = $array_ontology[0]->{'term_name'};
	
	$name =~ s/^\s+//;
	$name =~ s/\s+$//;
	
	$traits{ $trait } = $name;
	
#	print "$name -> $to -> $trait\n";
}

#my $phenotype_tree = '';
my $phenotype_tree = '<div class="label">Phenotypes:</div>'."\n";
$phenotype_tree .= '<div id="phenoTree" class="jstree"><ul>'."\n";
#$phenotype_tree .= '<li id="_at"><a href="#">All Traits</a><ul>'."\n";

my @keys = keys %traits;

for ( my $i = 0; $i < @keys; $i++ ) {
	$phenotype_tree .= "<li id=\"$keys[$i]\"><a href=\"#\">$keys[$i]: $traits{$keys[$i]}</a></li>\n";
}

#$phenotype_tree .= '</ul></li>'."\n";
$phenotype_tree .= '</ul></div>'."\n";

my($host,$user)=qw(127.0.0.1 root);
my($host,$user)=qw(brie.cshl.edu achuah);
my $self=$0;
$self=~s/^.+\///;
my @chrSize=(30427671,19698289,23459830,18585056,26975502);
my @cumSize=(0,0);
$cumSize[$_+1]=$cumSize[$_]+$chrSize[$_-1] for 1..$#chrSize;
unshift @chrSize,sum @chrSize;
my $chrSize=join',',@chrSize;
#print header('application/json');
print header('text/html');
my $q=new CGI;
autoEscape(0);
my $jsDir='';
$jsDir='/jshc' if $host=~/brie\.cshl/ || $host=~/dev\.gramene/;

my $dbh = Gramene::DB->new('diversity_arabidopsis');

#my($dsn, $password) = ("DBI:mysql:database=diversity_arabidopsis34;host=cabot;port=3306,;mysql_read_default_group=client");
#my $dbh = DBI->connect($dsn, $user, $password, {RaiseError=>1, AutoCommit=>1} ) or die "Could not connect.";

## Old dbh connect
#my $dbh=DBI->connect( "dbi:mysql:test;host=$host",$user,'',{ RaiseError => 1, AutoCommit => 0 } ) ||
#        die "Database connection not made: $DBI::errstr";

#my GWAS params
my $score=defined param('score') && param('score')?param('score'):'p_val'; # or rank or maf or mac
my $chr=defined param('chr') && param('chr')?param('chr'):0;
$chr=0 if $chr<0 || $chr>5;
my $plotLines='';
unless($chr) {
	$plotLines='plotLines: ['.join(",\n",map{"{color: '#000020', width: 2, value: $_}"} @cumSize[2..$#cumSize]).'],';
}
my $start=defined param('start') && param('start')?param('start')+0:1;
$start=1 if $start<1 || $start>$chrSize[$chr];
my $end=defined param('end') && param('end')?param('end')+0:$chrSize[$chr];
$end=$chrSize[$chr] if $end<1 || $end>$chrSize[$chr];
($start,$end)=($end,$start) if $start>$end;
my $rank=defined param('rank') && param('rank')?param('rank')+0:300;
$rank=1000 if $rank<1 || $rank>1000;
my $p_val=defined param('p_val') && param('p_val')?param('p_val'):0;
$p_val=0 if $p_val<0;
my $pheno=defined param('pheno') && param('pheno')?param('pheno'):'LD,SD,_fgx';
my $pheno0=$pheno;
$pheno=~s/_fgx/FRI,FLC/;
$pheno=~s/_mft/FT Diameter/;
$pheno=~s/_at1/At1,At1 CFU2/;
$pheno=~s/_at2/At2,At2 CFU2/;
$pheno=~s/_as2/As2,As2 CFU2/;
$pheno=~s/_as/As,As CFU2/;
$pheno=~s/_bs/Bs,Bs CFU2/;
$pheno=~s/_adr/Aphid/;
my @pheno=split',',$pheno;
my $phenoQ=join',',map{"'$_'"} @pheno;
my $open=defined param('open') && param('open')?param('open'):'_ft,_dtf,_fgx';
my $openQ=join',',map{"'$_'"} split',',$open;
my $axis='';
if($score eq 'p_val') {
        $axis='Score (-log10 p-value)';
}else {
        $axis=$score eq 'rank'?'Rank within Phenotype':$score eq 'maf'?'Minor Allele Frequency':'Minor Allele Coverage';
}

my $gwas_query = <<END;
    SELECT
        dtu.local_trait_name,
        cgs.chromosome,
        cgs.gwas_significance_value_blob,
        cgs.gwas_position_blob,
        cgs.gwas_rank_blob
    FROM
        cdv_g2p_study cgs,
        cdv_pheno_set cps,
        div_experiment de,
        div_obs_unit dou,
        div_trait dt,
        div_trait_uom dtu
    WHERE
        cgs.chromosome=?
    AND
        dtu.local_trait_name IN ($phenoQ)
    AND
        de.div_experiment_id=?
    AND
        de.div_experiment_id=dou.div_experiment_id
    AND
        dou.div_obs_unit_id=dt.div_obs_unit_id
    AND
        dt.div_trait_uom_id=dtu.div_trait_uom_id
    AND
        dt.div_trait_id=cps.div_trait_id
    AND
        cps.cdv_g2p_study_id=cgs.cdv_g2p_study_id
END

my $array_ref = $dbh->selectall_arrayref($gwas_query, { Columns => {} }, ( $chr, 5 ));

## Old select query
#my $where="and analysis='EmmaTrans' and gwas2010.call=32";
#$where.=" and chr='$chr'" if $chr;
#$where.=" and pos>$start" if $start>1;
#$where.=" and pos<$end" if $end<$chrSize[$chr];
#$where.=" and rank<=$rank" if $rank<1000;
#$where.=" and p_val>$p_val" if $p_val>2;
#my $sql=qq{select chr,pos,pheno,p_val,rank,maf,mac from gwas2010 where pheno in ($phenoQ) $where};
#my $results=$dbh->selectall_arrayref($sql,{Slice=>{}});
$dbh->disconnect();
my %x;

## Old data assignment
#for(@$results) {
#	push @{$x{$_->{'pheno'}}},[$_->{'chr'},$_->{'pos'},$_->{$score},$_->{'p_val'},$_->{'rank'},$_->{'maf'},$_->{'mac'}];
#}

for (@$array_ref) {
    my (
        $scores_header, @ranks
       ) = unpack "a200n*", $_->{'gwas_rank_blob'};
    my (
        $positions_header, @positions
       ) = unpack "a200N*", $_->{'gwas_position_blob'};
    my (
        $rank_header, @scores
       ) = unpack "a200f*", $_->{'gwas_significance_value_blob'};
    for ( my $i = 0; $i < @scores; $i++ ) {
        my $p_value = 10 ** -$scores[$i];
        if ((!$start or $positions[$i] >= $start) and (!$end or $positions[$i] <= $end) and (!$rank or $ranks[$i] <= $rank) and (!$p_val or $scores[$i] >= $p_val)) {
            push @{ $x{ $_->{'local_trait_name'} } }, [ $_->{'chromosome'}, $positions[$i], $scores[$i], $scores[$i], $ranks[$i], 10, 10 ];
        }
    }
}

my @series=();
my %order;
#my @color=();
#my @color=('223,83,83','119,152,191');
#my @color=('119,152,191','223,83,83','30,152,30','100,0,100'); #chloe's
my @color=('255,127,0','106,61,154','51,160,44','227,26,28','166,206,227','255,255,153','202,178,214','178,223,138','251,154,153','253,191,111','31,120,180','255,255,255');
my $i=0;
$order{$_}=++$i for @pheno;
$i=0;
for my $s(sort{$order{$a}<=>$order{$b}} keys %x) {
	my $dataString=join",\n",map{
			my($c,$x,$y,$p,$r,$maf,$mac)=@{$_};
			my $co='';
			unless($chr) {
				$co=",c:$c,co:$x";
				$x+=$cumSize[$c];
			};
			"{x:$x,y:$y,p:$p,r:$r,maf:$maf,mac:$mac$co}"
		} @{$x{$s}};
	my $defdColor=defined $color[$i]?",color:'rgba($color[$i],.6)'":'';
	push @series,qq({name:'$s'$defdColor,data:[$dataString]});
	$i++;
}
my $series=join(',',@series);
#die "$openQ $phenoQ" unless $pheno eq 'LD,SD,FRI,FLC';

print <<END;
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
	<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
	<title>Warelab's Arabidopsis Thaliana Association Study viz</title>
	<meta name="title" content="jsTree - 107 Atwell 2010 Phenotypes" /> 
	<script type="text/javascript" src="$jsDir/jsTree/_lib/jquery.js"></script>
<!--
	<script type="text/javascript" src="$jsDir/jquery/jquery-1.4.4.min.js"></script>
-->
	<link type="text/css" rel="stylesheet" href="$jsDir/jsTree/_docs/syntax/!style.css"/>
	<script type="text/javascript" src="$jsDir/jsTree/_docs/syntax/!script.js"></script>
	<link rel="stylesheet" type="text/css" href="$jsDir/jsTree/_docs/!style.css" />
	<script type="text/javascript" src="$jsDir/jsTree/_lib/jquery.cookie.js"></script>
	<script type="text/javascript" src="$jsDir/jsTree/_lib/jquery.hotkeys.js"></script>
	<script type="text/javascript" src="$jsDir/jsTree/jquery.jstree.js"></script>
	<script type="text/javascript" src="$jsDir/highcharts/js/highcharts.js"></script>
	<script type="text/javascript" src="$jsDir/highcharts/js/modules/exporting.js"></script>

<script type="text/javascript">
var example = 'scatter',
	theme = 'dark-blue';
SyntaxHighlighter.config.clipboardSwf = "$jsDir/jsTree/_docs/syntax/clipboard.swf";	
</script>
	
	<link href="$jsDir/jquery-ui/themes/base/jquery.ui.all.css" rel="stylesheet" type="text/css"/>
	<script type="text/javascript" src="$jsDir/splitter/splitter.js"></script>
	<style type="text/css" media="all">
html, body{
height:100%;width:100%; 
margin:0;padding:0;overflow: hidden;
}
#header{background:#c4dcfb;height:26px;}/* Header */

#splitterContainer {/* Main splitter element */
height:97%;width:100%;margin:0;border:0;
}

#leftPane{
float:left;width:75%;height:100%;border-top:solid 1px #9cbdff;
background:#000000;
}
#rightPane{	/*Contains toolbar and horizontal splitter*/
float:right;width:25%;height:100%;
background:#FFFFEE;
overflow: auto;
}

/* Splitbar styles; these are the default class names and required styles */
.splitbarV {
float:left;width:6px;height:100%;
line-height:0px;font-size:0px;
border-left:solid 1px #9cbdff;border-right:solid 1px #9cbdff;
background:#cbe1fb url($jsDir/splitter/img/panev.gif) 0% 50%;
}
.splitbarH {
height:6px;text-align:left;line-height:0px;font-size:0px;
border-top:solid 1px #9cbdff;border-bottom:solid 1px #9cbdff;
background:#cbe1fb url($jsDir/splitter/img/paneh.gif) 50% 0%;
}

.splitbuttonV{
margin-top:-41px;margin-left:-4px;top:50%;position:relative;
height:83px;width:10px;
background:transparent url($jsDir/splitter/img/panevc.gif) 10px 50%;
}
.splitbuttonV.invert{
margin-left:0px;background:transparent url($jsDir/splitter/img/panevc.gif) 0px 50%;
}
.splitbuttonH{
margin-left:-41px;left:50%;position:relative;
height:10px !important;width:83px;
background:transparent url($jsDir/splitter/img/panehc.gif) 50% 0px;
}
.splitbuttonH.invert{
margin-top:-4px;background:transparent url($jsDir/splitter/img/panehc.gif) 50% -10px;
}
.splitbarV.working,.splitbarH.working,.splitbuttonV.working,.splitbuttonH.working{
 -moz-opacity:.50; filter:alpha(opacity=50); opacity:.50;
}

#demo-frame > div.header { padding: 10px !important; }
</style>
	
<script type="text/javascript" class="source">

var chrSize=[$chrSize];
var chr=$chr;
var chrLabel='Chr '+chr;
if(chr==0) { chrLabel='Genome'; }
var label={ 'p_val':'p', 'rank':'Rank', 'maf':'maf', 'mac':'mac' };
var chart;

var options= {
	colors: ["#DDDF0D", "#55BF3B", "#DF5353", "#7798BF", "#aaeeee", "#ff0066", "#eeaaee", "#55BF3B", "#DF5353", "#7798BF", "#aaeeee"],
	chart: {
		renderTo: 'container',
		defaultSeriesType: 'scatter',
		zoomType: 'xy',
		backgroundColor: {
			linearGradient: [0, 0, 250, 500],
			stops: [
				[0,'rgb(48, 48, 96)'],
				[1,'rgb(0, 0, 0)']
			]
		},
		borderColor: '#000000',
		className: 'dark-container',
		plotBackgroundColor: 'rgba(255, 255, 255, .1)',
		plotBorderColor: '#CCCCCC',
		plotShadow: false
	},
	title: {
		text: 'Arabidopsis thaliana 107 Phenotype Association',
		style: {
			color: '#C0C0C0'
		}
	},
	subtitle: {
		text: 'Source: Atwell et.al., Nature 2010'
	},
	xAxis: {
		title: {
			enabled: true,
			text: chrLabel+' pos (Mb)',
			style: {
				color: '#C0C0C0'
            }
		},
		labels: {
            rotation: -45,
            align: 'right',
			backgroundColor: '#FFFFFF',
			style: {
				color: '#A0A0A0'
			},
            formatter: function() {
                return this.value/1e6;
            }
        },
        $plotLines
		gridLineWidth: 1,
		gridLineColor: '#333333',
		lineColor: '#A0A0A0',
		tickColor: '#A0A0A0'
	},
	yAxis: {
		title: {
			text: '$axis',
			style: {
				color: '#C0C0C0'
			}
		},
		gridLineColor: '#333333',
		labels: {
			style: {
				color: '#A0A0A0'
			}
		},
		lineColor: '#A0A0A0',
		minorTickInterval: null,
		tickColor: '#A0A0A0'
	},
	tooltip: {
		formatter: function() {
			var num=Math.pow(10,-this.y);
			var c=chr;
			var co=this.x;
			if(c<1) {
				c=this.point.c;
				co=this.point.co;
			}
			return 'Chr '+c+':'+co+'<br> p = '+Math.pow(10,-this.point.p).toPrecision(4)+'<br> maf = '+this.point.maf.toPrecision(4)+'<br> mac = '+this.point.mac+'<br>'+this.series.name+' #'+this.point.r;
		},
		backgroundColor: 'rgba(0, 0, 0, 0.75)',
		style: {
			color: '#F0F0F0'
		}
	},
	toolbar: {
		itemStyle: { 
			color: 'silver'
		}
	},	
	legend: {
		floating: true,
		layout: 'vertical',
		align: 'right',
		verticalAlign: 'top',
		x: -3,
		y: 3,
		backgroundColor: 'rgba(0, 0, 0, 0.5)',
		borderWidth: 1,
		style: {
			color: '#A0A0A0'
		},
		itemStyle: {
			color: '#CCC'
		},
		itemHoverStyle: {
			color: '#FFF'
		},
		itemHiddenStyle: {
			color: '#444'
		}
	},
	credits: {
		enabled: false,
		style: {
			color: '#666'
		}
	},
	labels: {
		style: {
			color: '#CCC'
		}
	},
	legendBackgroundColor: 'rgba(0, 0, 0, 0.5)',
	legendBackgroundColorSolid: 'rgb(35, 35, 70)',
	dataLabelsColor: '#444',
	maskColor: 'rgba(255,255,255,0.3)',
	exporting: {
		enabled: true,
		buttons: {
			exportButton: {
				verticalAlign: 'bottom',
				y: -4
			},
			printButton: {
				enabled: false,
				verticalAlign: 'bottom',
				x: -33,
				y: -5
			}
		}
	},
	plotOptions: {
		scatter: {
			marker: {
				radius: 5,
				states: {
					hover: {
						enabled: true,
						lineColor: 'rgb(100,100,100)'
					}
				}
			},
			states: {
				hover: {
					marker: {
						enabled: false
					}
				}
			}
		},					
		series: {
			point: {
				events: {
					click: function() {
						var c=chr;
						var co=this.x;
						if(c<1) {
							c=this.c;
							co=this.co;
						}
						if(this.series.data.length>1) {
							window.open("http://www.gramene.org/Arabidopsis_thaliana/Location/View?r="+c+":"+(co-100)+"-"+(co+99));
						}
					}
				}
			}
		}
	},
	series: [$series]
};
/*

*/



\$(document).ready(function () {
	chart=new Highcharts.Chart(options);
	\$("#splitterContainer").splitter({minBsize:0,maxBsize:400,splitVertical:false,A:\$('#leftPane'),B:\$('#rightPane'),closeableto:100});
	\$("#phenoTree").jstree({
		"plugins" : [ "themes", "html_data", "cookies", "checkbox", "ui" ],
		"cookies" : { "save_selected": false },
		"ui" : { "initially_select" : [ $phenoQ ] },
		"core" : { "initially_open" : [ $openQ ] },
		"callback" : {
			"onload" : function (tree) {
				\$('li[selected=true]').each(function () {
					\$.tree.plugins.checkbox.check(this);
				});
			}
		}
	});
	\$("#phenoTree").bind("select_node.jstree", function(e, data) {
		data.rslt.obj.parents('.jstree-closed').each(function() {
			data.inst.open_node(this);
  		});
	});
	\$("#phenoTree").click(function() {
		var checked_ids = [];
		\$('#phenoTree').jstree('get_checked').each(function () {
			checked_ids.push(this.id);
		});
/*
		jQuery.jstree._focused().get_checked().each(function () {
			checked_ids.push(this.id);
		});
		var selecookie=\$.cookie("jstree_select");
		var opencookie=\$.cookie("jstree_open");
		selecookie=selecookie.replace(/#/g,'');
		opencookie=opencookie.replace(/#/g,'');
		\$("#open").val(opencookie);
*/		
		var selected=checked_ids.join(",");
		if(selected!="$pheno0") {
			\$("#pheno").val(selected);
			\$('#form').submit();
		}
	});
	\$("#chr").change(function() { 
        \$("#start").val(1);
        \$("#end").val(chrSize[\$("#chr").val()]);
		\$('#form').submit();
	});
	\$("#score").change(function() {
		\$('#form').submit();
	});
});

// Start Automatic resize mod
		var containerWidth,containerHeight;			
		function checkContainerResize() {
			var \$container = \$('#container'),
				width = \$container.width(),
				height = \$container.height();
			if (containerWidth === undefined || containerHeight == undefined) {
				containerWidth = width;
				containerHeight = height;
			} 
			if (width != containerWidth || height != containerHeight) {				
				containerWidth = width;
				containerHeight = height;
				reflowChart();
			}
		}		
		function reflowChart () {
			var i = options.series.length;
			while (i--) {
				options.series[i].animation = false; 
			}
			chart=new Highcharts.Chart(options);
		}
		\$(window).load(function() {
			setInterval(checkContainerResize, 1000);
		});
// End automatic resize mod 
</script>	
	
</head>
<body>
<div id="header">
	<form id="form" method="post" action="$self">
	<input id="pheno" name="pheno" type="hidden" value="$pheno">
	<input id="open" name="open" type="hidden" value="$open">
	<select id="score" name="score">
		<option value="p_val"@{[$score eq "p_val"?' selected':'']}>p-value</option>
		<option value="rank"@{[$score eq "rank"?' selected':'']}>Rank</option>
		<option value="maf"@{[$score eq "maf"?' selected':'']}>MAF</option>
		<option value="mac"@{[$score eq "mac"?' selected':'']}>MAC</option>
	</select>
	<select id="chr" name="chr"> 
<!--		<option value="0"@{[$chr==0?' selected':'']}>Genome</option> -->
		<option value="1"@{[$chr==1?' selected':'']}>Chr 1</option>
		<option value="2"@{[$chr==2?' selected':'']}>Chr 2</option>
		<option value="3"@{[$chr==3?' selected':'']}>Chr 3</option>
		<option value="4"@{[$chr==4?' selected':'']}>Chr 4</option> 
		<option value="5"@{[$chr==5?' selected':'']}>Chr 5</option>
	</select>:<input id="start" name="start" type="text" value="$start" maxlength=11 size=6>-<input id="end" name="end" type="text" value="$end" maxlength=11 size=6>
	Filter:&nbsp;Top<input id="rank" name="rank" value=$rank maxlength=4 size=2>SNPs/pheno, 
	Score&ge;<input id="p_val" name="p_val" value="$p_val" maxlength=4 size=2>
	<input type="submit" value="Go">
	</form>
</div>
<div id="splitterContainer">
	<div id="leftPane">
		<div id="container" style="width: 100%; height: 100%; margin: 0 auto"></div>
	</div>
	<div id="rightPane" style="margin: 0 auto">
$phenotype_tree
</div>
</body>
</html>

END
1;

