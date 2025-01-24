#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;

#File: run_container.pl
#Description: Runs container modeling program for specified experiments.  Constructs parameter file based on a template, then runs program.  Runs plots.
#Input: Container program, specify parameters here, input met data
#Output: Container program output files
#Dan Steinhoff
# 13 March 2013

# Revision history (CNY)
#
# changed the following user specifications:
# (make sure to change these when changing cities AND when changing timespans!!!!)
#   location variables
#   manual fill flag: off
#   whs_mul (initial water height multiplier): set to 0.1 (originally 0.0)
#   tws (initial water temp): set to the initial air temp (T) value in each city's text file
#   tgs (initial ground temp): set to the initial soil temp (TSOIL) value in each city's text file
#   cflag (whether clouds specified): set to "F" when using GEOS S2S input
#   sflag (whether soil temp specified): set to "F" when using GEOS S2S input
# added more temperature bins for frequency stats
# added a run name that summarizes a few specs: city, timespan, manual fill (this is also used for generating directory name)
# changed ncl plotting to remove ncl copyright info from showing up every time (changed "ncl ..." to "ncl -Q ...")
# 2023-08-09, CNY: added manual emptying

#### USER SPECIFICATIONS ####
# Run name
my $dataSrc = "GEOS-S2S_"; # "";
my $runName = "${dataSrc}Ng_Ja_NE_2017100100_2017103123_noFill"; # Summary of specifications
# Location
my @cities = qw/ Negombo Jaffna NuwaraEliya /; # List of cities ; Trincomalee
my @sublocs = qw/ 2017100100_2017103123 2017100100_2017103123 2017100100_2017103123 /; # List of timespans for cities (YYYYMMDDHH) (a bit of a hack: using the sublocation string as the timespan) (MUST match input txt filenames)
my @lats  = (7.2008,  9.6615,  6.9497);   # Site latitudes (deg N)  ; 8.5874 for Trinco
my @lons  = (79.8737, 80.0255, 80.7891);  # Site longitudes (deg E) ; 81.2152 for Trinco
my @elevs = (2.,      5.,      1868.);    # Site elevations (m)     ; 8. for Trinco
# Containers
my @cshapes  = qw/ ROUND /;    # qw/ ROUND ROUND ROUND /;  # Shape of containers ('ROUND', 'RECT', 'TIRE')
my @radii1_t = (0.131);        # (0.082, 0.131, 0.298);    # Top container lengths 1 (in m)
my @radii2_t = (0.131);        # (0.082, 0.131, 0.298);    # Top container lengths 2 (in m)
my @radii1_b = (0.131);        # (0.082, 0.131, 0.298);    # Body container lengths 1 (in m)
my @radii2_b = (0.131);        # (0.082, 0.131, 0.298);    # Body container lengths 2 (in m)
my @heights  = (0.368);        # (0.191, 0.368, 0.921);    # Container heights (in m)
my @conducs  = (0.50);         # (0.30,  0.50,  0.50);     # Thermal conductivity of containers (in W m-1 K-1)
my @thicks   = (0.0023);       # (0.0018, 0.0023, 0.0048); # Thickness of containers (in m)
my @names_cont = qw/ Bucket /; # qw/ Coffee Bucket Barrel /; # Names of containers
# Reflectivities (colors)
my @refls = (0.5);             # (0.1, 0.5, 0.9); # Reflectivity coefficients of containers - depends on number of containers
my @names_refl = qw/ Gray /;   # qw/ Black Gray White /; # Names of reflectivities (colors) - depends on number of containers
# Shade
my @shades = (0.50);           # (0.00, 0.50, 1.00); # Fraction of shade
my @names_shade = qw/ HalfShade /; # qw/ NoShade HalfShade FullShade /; # Names of shade conditions
# Clouds (if not provided)
my @clouds = ([0.0,0.0,0.0,0.0,0.0,0.0], [0.33,0.33,0.33,0.1,0.1,0.1], [0.66,0.66,0.66,0.33,0.33,0.33]); # Inner index is low, middle, and high clouds, for day and night.
my @names_cloud = qw/ Clear PartlyCloudy Overcast /; # Names of cloud cover conditions (if specified)
# Water, ground, and cloud conditions
my @whs_mul = (0.1, 0.1, 0.1); # Multiplier of initial water height / container height ratio - depends on number of containers
my @tws = (25.37, 27.90, 20.75); # Initial water temperatures (deg C)  - depends on number of locations (for 2014-15: Ng 23.88, Ja 26.66, NE 19.57, Tr 25.89; for 2016-17: Ng 23.42, Ja 24.84, NE 16.71, Tr 24.48)
my @tgs = (26.30, 30.65, 22.67); # Initial ground temperatures (deg C) - depends on number of locations (for 2014-15: Ng 22.73, Ja 25.11, NE 19.87, Tr 24.45; for 2016-17: Ng 22.17, Ja 23.94, NE 18.49, Tr 23.78)
# Flags
my $cflag = "F"; # Whether cloud fraction data is input (T/F)    !! set to "T" when using MERRA/IMERG, "F" when using GEOS S2S
my $sflag = "F"; # Whether soil temperature data is input (T/F)  !! set to "T" when using MERRA/IMERG, "F" when using GEOS S2S
my $mfill = "F"; # Manual fill flag (T/F)
my $fill_intv = 168; # Manual fill interval, in hours
my $mempty = "F"; # Manual empty flag (T/F)               (ADDED ON 2023-08-09, CNY)
my $empty_intv = 168; # Manual empty interval, in hours   (ADDED ON 2023-08-09, CNY)

# Misc.
my $time_step = 60.00; # Model time step in seconds
my $input_intv = 3600.00; # Time interval of input data in seconds
my $miss = -9999.00; # Missing value of input data
my $plot = 0; # Whether to run plotting routines or not (1/0)  # set this to zero if doing only one container, as plotting expects data for each container type
my $path = "/home/local/WIN/cyasana1/Code/whatchem-2.1_CNY/container_modeling/src"; # Path to container.exe file  ***
my $path_out_data = "/mnt/redwood/local_drive/chanud/Output/whatchem_2.1_CNY/data"; # Path for output data
my $path_out_plots = "/mnt/redwood/local_drive/chanud/Output/whatchem_2.1_CNY/plots"; # Path for output plots
#### END USER SPECIFICATIONS ####

# Welcome message
print "##############################\n";
print "#### WELCOME TO WHATCH'EM ####\n";
print "##############################\n";
print "# Starting...\n";
print "# Run name: ${runName}\n";

# Get array sizes for use throughout program
my $exp_cont_size = scalar (@radii1_t);
my $exp_shade_size = scalar (@shades);
my $exp_refl_size = scalar (@refls);
my $exp_cloud_size = scalar (@clouds);

# Open master statistics files
open MTOUTFILE, ">stats/stats_temperature.txt" or die "Cannot open master temperature stats file: $!";
open MWOUTFILE, ">stats/stats_waterheight.txt" or die "Cannot open master water height stats file: $!";
if ($cflag eq 'F')
{
	printf MTOUTFILE "%-30s %-10s %-15s %-10s %-15s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s\n", "CITY-HOBO", "CONTAINER", "REFLS", "SHADE", "CLOUDS", "TMEAN", "TMAX", "TMIN", "CMEAN", "CMAX", "CMIN", "TAMEAN", "EVAP";
	printf MWOUTFILE "%-30s %-10s %-15s %-10s %-15s %-15s\n", "CITY-HOBO", "CONTAINER", "REFLS", "SHADE", "CLOUDS", "WATER HEIGHT";
}
else
{
	printf MTOUTFILE "%-30s %-10s %-15s %-10s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s\n", "CITY-HOBO", "CONTAINER", "REFLS", "SHADE", "TMEAN", "TMAX", "TMIN", "CMEAN", "CMAX", "CMIN", "TAMEAN", "EVAP";
	printf MWOUTFILE "%-30s %-10s %-15s %-10s %-15s\n", "CITY-HOBO", "CONTAINER", "REFLS", "SHADE", "WATER HEIGHT";
}

# Loop through each city
my @sublocs_t = @sublocs;
foreach my $city (@cities)
{
	print "# $city\n";
	# Get the city-specific variables
	my $subloc = shift (@sublocs_t);
	my $input_file = "input/container_input_${dataSrc}$city-$subloc.txt";
	print "# Input file: $input_file\n";
	my $output_file_base = "container_output_${dataSrc}$city-$subloc";
	if (!(-d "$path_out_data/${dataSrc}$city-$subloc")) { system "mkdir $path_out_data/${dataSrc}$city-$subloc" }
	if (!(-d "$path_out_plots/${runName}")) # new
	{ 
		system "mkdir $path_out_plots/${runName}";
		system "mkdir $path_out_plots/${runName}/avg_temperature";
		system "mkdir $path_out_plots/${runName}/bar_dtr";
		system "mkdir $path_out_plots/${runName}/bar_temperature";
		system "mkdir $path_out_plots/${runName}/bar_waterheight";
		system "mkdir $path_out_plots/${runName}/clouds_precip";
		system "mkdir $path_out_plots/${runName}/comparison_cities";
		system "mkdir $path_out_plots/${runName}/comparison_clouds";
		system "mkdir $path_out_plots/${runName}/comparison_containers";
		system "mkdir $path_out_plots/${runName}/comparison_refls";
		system "mkdir $path_out_plots/${runName}/comparison_shade";
		system "mkdir $path_out_plots/${runName}/energy_balance";
		system "mkdir $path_out_plots/${runName}/freq_distribution";
		system "mkdir $path_out_plots/${runName}/temp_and_rh";
		system "mkdir $path_out_plots/${runName}/temperature";
	}
	my $lat = sprintf "%.2f", shift (@lats);
	my $lon = sprintf "%.2f", shift (@lons);
	my $elev = sprintf "%.2f", shift (@elevs);
	my $tw = sprintf "%.2f", shift (@tws);
	my $tg = sprintf "%.2f", shift (@tgs);
	my $t_diff = $input_intv/$time_step;
	$time_step = sprintf "%.2f", $time_step;
	$t_diff = sprintf "%.2f", $t_diff;

	# Loop through each container
	for (my $i = 0; $i < $exp_cont_size; $i++)
	{
		my $cshape = $cshapes[$i];
		my $radius1_t = sprintf "%.4f", $radii1_t[$i];
		my $radius2_t = sprintf "%.4f", $radii2_t[$i];
		my $radius1_b = sprintf "%.4f", $radii1_b[$i];
		my $radius2_b = sprintf "%.4f", $radii2_b[$i];
		my $height = sprintf "%.4f", $heights[$i];
		my $wh = sprintf "%.4f", $height*$whs_mul[$i];
		my $name_cont = $names_cont[$i];
		print "## $name_cont\n";

		# Loop through each reflectivity
		for (my $j = 0; $j < $exp_refl_size; $j++)
		{
			my $name_refl = $names_refl[$j];
			print "### $name_refl\n";
			my $refl = sprintf "%.2f", $refls[$j];
			my $conduc = sprintf "%.2f", $conducs[$j];
			my $thick = sprintf "%.4f", $thicks[$j];

			# Loop through each shade
			for (my $k = 0; $k < $exp_shade_size; $k++)
			{
				my $name_shade = $names_shade[$k];
				print "#### $name_shade\n";
				my $shade = sprintf "%.2f", $shades[$k];

				# Loop through each cloud (if not provided)
				if ($cflag eq 'F')
				{
					for (my $l = 0; $l < $exp_cloud_size; $l++)
					{
						my $name_cloud = $names_cloud[$l];
						print "##### $name_cloud\n";
						my $output_file = $output_file_base."-$name_cont-$name_refl-$name_shade-$name_cloud.txt";
						my @clouds_exp = @{$clouds[$l]};
						my $cld = sprintf "%.2f", $clouds_exp[0];
						my $cmd = sprintf "%.2f", $clouds_exp[1];
						my $chd = sprintf "%.2f", $clouds_exp[2];
						my $cln = sprintf "%.2f", $clouds_exp[3];
						my $cmn = sprintf "%.2f", $clouds_exp[4];
						my $chn = sprintf "%.2f", $clouds_exp[5];

						# Construct parameter file
						print "##### Constructing parameter file...\n";
						open OUTFILE, ">param_file.txt";
						print OUTFILE "$input_file\n";
						print OUTFILE "$city\n";
						print OUTFILE "$subloc\n";
						print OUTFILE "$lat\n";
						print OUTFILE "$lon\n";
						print OUTFILE "$elev\n";
						print OUTFILE "$cshape\n";
						print OUTFILE "$radius1_t\n";
						print OUTFILE "$radius2_t\n";
						print OUTFILE "$radius1_b\n";
						print OUTFILE "$radius2_b\n";
						print OUTFILE "$height\n";
						print OUTFILE "$refl\n";
						print OUTFILE "$conduc\n";
						print OUTFILE "$thick\n";
						print OUTFILE "$wh\n";
						print OUTFILE "$tw\n";
						print OUTFILE "$tg\n";
						print OUTFILE "$time_step\n";
						print OUTFILE "$t_diff\n";
						print OUTFILE "$shade\n";
						print OUTFILE "$cflag\n";
						print OUTFILE "$sflag\n";
						print OUTFILE "$mfill\n";
						print OUTFILE "$fill_intv\n";
						print OUTFILE "$mempty\n";       # (ADDED ON 2023-08-09, CNY)
						print OUTFILE "$empty_intv\n";   # (ADDED ON 2023-08-09, CNY)
						print OUTFILE "$output_file\n";
						print OUTFILE "$cld\n";
						print OUTFILE "$cmd\n";
						print OUTFILE "$chd\n";
						print OUTFILE "$cln\n";
						print OUTFILE "$cmn\n";
						print OUTFILE "$chn\n";
						close OUTFILE;

						# Run container modeling program
						print "##### Running container model...";
						system "$path/container.exe";

						# Move output file to correct directory
						print "##### Writing output files...\n";
						system "mv $output_file $path_out_data/${dataSrc}$city-$subloc";

						# Run daily energy balance and basic stats
						my $einfile = "$path_out_data/${dataSrc}$city-$subloc/$output_file";
						my $eoutfile = "$path_out_data/${dataSrc}$city-$subloc/energybal_avg_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade-$name_cloud.txt";
						my $soutfile = "$path_out_data/${dataSrc}$city-$subloc/stats_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade-$name_cloud.txt";
						&daily_energy_balance_stats($einfile, $eoutfile, $soutfile, "$city-$subloc", $name_cont, $name_refl, $name_shade, $cflag, $name_cloud);

						# Run first set of plots
						if ($plot)
						{
							print "##### Running temperature and energy balance plots...\n";
							# Plotting - Temperature variables (tw, ta, tg)
							system "ncl -Q \'infile=\"$path_out_data/${dataSrc}$city-$subloc/$output_file\"\' \'outfile=\"$path_out_plots/${runName}/temperature/$city-$subloc-$name_cont-$name_refl-$name_shade-$name_cloud\"\' \'container=\"$name_cont\"\' \'refl=\"$name_refl\"\' \'shade=\"$name_shade\"\' \'clouds=\"$name_cloud\"\' \'site=\"$city-$subloc\"\' plots/plot_t.ncl";

							# Plotting - Average energy balance - Water
							system "ncl -Q \'infile=\"$path_out_data/${dataSrc}$city-$subloc/energybal_avg_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade-$name_cloud.txt\"\' \'outfile1=\"$path_out_plots/${runName}/energy_balance/$city-$subloc-$name_cont-$name_refl-$name_shade-water\"\' \'outfile2=\"$path_out_plots/${runName}/energy_balance/$city-$subloc-$name_cont-$name_refl-${name_shade}\"\' \'container=\"$name_cont\"\' \'refl=\"$name_refl\"\' \'shade=\"$name_shade\"\' \'clouds=\"$name_cloud\"\' \'site=\"$city-$subloc\"\' plots/plot_energy_balance_avg_water.ncl";

							# Plotting - Average energy balance - Container
                            system "ncl -Q \'infile=\"$path_out_data/${dataSrc}$city-$subloc/energybal_avg_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade-$name_cloud.txt\"\' \'outfile=\"$path_out_plots/${runName}/energy_balance/$city-$subloc-$name_cont-$name_refl-$name_shade-container\"\' \'container=\"$name_cont\"\' \'refl=\"$name_refl\"\' \'shade=\"$name_shade\"\' \'clouds=\"$name_cloud\"\' \'site=\"$city-$subloc\"\' plots/plot_energy_balance_avg_con.ncl";

							# Plotting = Average daily water temperature
							system "ncl -Q \'infile=\"$path_out_data/${dataSrc}$city-$subloc/energybal_avg_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade-$name_cloud.txt\"\' \'outfile=\"$path_out_plots/${runName}/avg_temperature/$city-$subloc-$name_cont-$name_refl-$name_shade-avg\"\' \'container=\"$name_cont\"\' \'refl=\"$name_refl\"\' \'shade=\"$name_shade\"\' \'clouds=\"$name_cloud\"\' \'site=\"$city-$subloc\"\' plots/plot_t_avg.ncl";
						}
					}
					# Plotting - Cloud experiments
					if ($plot)
					{
						print "##### Running cloud experiment plots...\n";
						system "ncl -Q \'infile1=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade-Clear.txt\"\' \'infile2=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade-PartlyCloudy.txt\"\' \'infile3=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade-Overcast.txt\"\' \'outfile1=\"$path_out_plots/${runName}/comparison_clouds/temperature-$city-$subloc-$name_cont-$name_refl-$name_shade-clouds\"\' \'outfile2=\"$path_out_plots/${runName}/comparison_clouds/waterheight-$city-$subloc-$name_cont-$name_refl-$name_shade-clouds\"\' \'container=\"$name_cont\"\' \'refl=\"$name_refl\"\' \'shade=\"$name_shade\"\' \'site=\"$city-$subloc\"\' plots/plot_comparison_cloud.ncl";
					}
				}
				else
				{
					my $output_file = $output_file_base."-$name_cont-$name_refl-$name_shade.txt";
					# Construct parameter file
					print "#### Constructing parameter file...\n";
					open OUTFILE, ">param_file.txt";
					print OUTFILE "$input_file\n";
					print OUTFILE "$city\n";
					print OUTFILE "$subloc\n";
					print OUTFILE "$lat\n";
					print OUTFILE "$lon\n";
					print OUTFILE "$elev\n";
					print OUTFILE "$cshape\n";
					print OUTFILE "$radius1_t\n";
					print OUTFILE "$radius2_t\n";
					print OUTFILE "$radius1_b\n";
					print OUTFILE "$radius2_b\n";
					print OUTFILE "$height\n";
					print OUTFILE "$refl\n";
					print OUTFILE "$conduc\n";
					print OUTFILE "$thick\n";
					print OUTFILE "$wh\n";
					print OUTFILE "$tw\n";
					print OUTFILE "$tg\n";
					print OUTFILE "$time_step\n";
					print OUTFILE "$t_diff\n";
					print OUTFILE "$shade\n";
					print OUTFILE "$cflag\n";
					print OUTFILE "$sflag\n";
					print OUTFILE "$mfill\n";
					print OUTFILE "$fill_intv\n";
					print OUTFILE "$mempty\n";       # (ADDED ON 2023-08-09, CNY)
					print OUTFILE "$empty_intv\n";   # (ADDED ON 2023-08-09, CNY)
					print OUTFILE "$output_file\n";
					close OUTFILE;

					# Run container modeling program
					print "#### Running container model...";
					system "$path/container.exe";

					# Move output file to correct directory
					print "#### Writing output files...\n";
					system "mv $output_file $path_out_data/${dataSrc}$city-$subloc";

					# Run daily energy balance and basic stats
					my $einfile = "$path_out_data/${dataSrc}$city-$subloc/$output_file";
					my $eoutfile = "$path_out_data/${dataSrc}$city-$subloc/energybal_avg_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade.txt";
					my $soutfile = "$path_out_data/${dataSrc}$city-$subloc/stats_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade.txt";
					&daily_energy_balance_stats($einfile, $eoutfile, $soutfile, "$city-$subloc", $name_cont, $name_refl, $name_shade, $cflag);

					# Plotting routines
					if ($plot)
					{
						print "#### Running temperature and energy balance plots...\n";
						# Plotting - Temperature variables (tw, ta, tg)
						system "ncl -Q \'infile=\"$path_out_data/${dataSrc}$city-$subloc/$output_file\"\' \'outfile=\"$path_out_plots/${runName}/temperature/$city-$subloc-$name_cont-$name_refl-$name_shade\"\' \'container=\"$name_cont\"\' \'refl=\"$name_refl\"\' \'shade=\"$name_shade\"\' \'clouds=\"Observed\"\' \'site=\"$city-$subloc\"\' plots/plot_t.ncl";

						# Plotting - Average energy balance - Water
						system "ncl -Q \'infile=\"$path_out_data/${dataSrc}$city-$subloc/energybal_avg_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade.txt\"\' \'outfile1=\"$path_out_plots/${runName}/energy_balance/$city-$subloc-$name_cont-$name_refl-$name_shade-water\"\' \'outfile2=\"$path_out_plots/${runName}/energy_balance/$city-$subloc-$name_cont-$name_refl-$name_shade\"\' \'container=\"$name_cont\"\' \'refl=\"$name_refl\"\' \'shade=\"$name_shade\"\' \'clouds=\"Observed\"\' \'site=\"$city-$subloc\"\' plots/plot_energy_balance_avg_water.ncl";

						# Plotting - Average energy balance - Container
                        system "ncl -Q \'infile=\"$path_out_data/${dataSrc}$city-$subloc/energybal_avg_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade.txt\"\' \'outfile=\"$path_out_plots/${runName}/energy_balance/$city-$subloc-$name_cont-$name_refl-$name_shade-container\"\' \'container=\"$name_cont\"\' \'refl=\"$name_refl\"\' \'shade=\"$name_shade\"\' \'clouds=\"Observed\"\' \'site=\"$city-$subloc\"\' plots/plot_energy_balance_avg_con.ncl";

						# Plotting = Average daily water temperature
						system "ncl -Q \'infile=\"$path_out_data/${dataSrc}$city-$subloc/energybal_avg_${dataSrc}$city-$subloc-$name_cont-$name_refl-$name_shade.txt\"\' \'outfile=\"$path_out_plots/${runName}/avg_temperature/$city-$subloc-$name_cont-$name_refl-$name_shade-avg\"\' \'container=\"$name_cont\"\' \'refl=\"$name_refl\"\' \'shade=\"$name_shade\"\' \'clouds=\"Observed\"\' \'site=\"$city-$subloc\"\' plots/plot_t_avg.ncl";
					}
				}
			}
		}

		if ($plot)
		{
			# Plotting - Reflectivity experiments
			print "#### Running color experiment plots...\n";
			for (my $k = 0; $k < $exp_shade_size; $k++)
			{
				my $name_shade = $names_shade[$k];
				if ($cflag eq 'F')
				{
					for (my $l = 0; $l < $exp_cloud_size; $l++)
					{
						my $name_cloud = $names_cloud[$l];
						system "ncl -Q \'infile1=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-Black-$name_shade-$name_cloud.txt\"\' \'infile2=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-Gray-$name_shade-$name_cloud.txt\"\' \'infile3=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-White-$name_shade-$name_cloud.txt\"\' \'outfile1=\"$path_out_plots/${runName}/comparison_refls/temperature-$city-$subloc-$name_cont-$name_shade-$name_cloud-refls\"\' \'outfile2=\"$path_out_plots/${runName}/comparison_refls/waterheight-$city-$subloc-$name_cont-$name_shade-$name_cloud-refls\"\' \'container=\"$name_cont\"\' \'shade=\"$name_shade\"\' \'clouds=\"$name_cloud\"\' \'site=\"$city-$subloc\"\' plots/plot_comparison_refl.ncl";
					}
				}
				else
				{
					system "ncl -Q \'infile1=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-Black-$name_shade.txt\"\' \'infile2=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-Gray-$name_shade.txt\"\' \'infile3=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-White-$name_shade.txt\"\' \'outfile1=\"$path_out_plots/${runName}/comparison_refls/temperature-$city-$subloc-$name_cont-$name_shade-refls\"\' \'outfile2=\"$path_out_plots/${runName}/comparison_refls/waterheight-$city-$subloc-$name_cont-$name_shade-refls\"\' \'container=\"$name_cont\"\' \'shade=\"$name_shade\"\' \'clouds=\"Observed\"\' \'site=\"$city-$subloc\"\' plots/plot_comparison_refl.ncl";
				}
			}
			# Plotting - Shade experiments
			print "### Running shade experiment plots...\n";
			for (my $j = 0; $j < $exp_refl_size; $j++)
			{
				my $name_refl = $names_refl[$j];
				if ($cflag eq 'F')
				{
					for (my $l = 0; $l < $exp_cloud_size; $l++)
					{
						my $name_cloud = $names_cloud[$l];
						system "ncl -Q \'infile1=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-$name_refl-NoShade-$name_cloud.txt\"\' \'infile2=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-$name_refl-HalfShade-$name_cloud.txt\"\' \'infile3=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-$name_refl-FullShade-$name_cloud.txt\"\' \'outfile1=\"$path_out_plots/${runName}/comparison_shade/temperature-$city-$subloc-$name_cont-$name_refl-$name_cloud-shade\"\' \'outfile2=\"$path_out_plots/${runName}/comparison_shade/waterheight-$city-$subloc-$name_cont-$name_refl-$name_cloud-shade\"\' \'container=\"$name_cont\"\' \'refl=\"$name_refl\"\' \'clouds=\"$name_cloud\"\' \'site=\"$city-$subloc\"\' plots/plot_comparison_shade.ncl";
					}
				}
				else
				{
					system "ncl -Q \'infile1=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-$name_refl-NoShade.txt\"\' \'infile2=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-$name_refl-HalfShade.txt\"\' \'infile3=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-$name_cont-$name_refl-FullShade.txt\"\' \'outfile1=\"$path_out_plots/${runName}/comparison_shade/temperature-$city-$subloc-$name_cont-$name_refl-shade\"\' \'outfile2=\"$path_out_plots/${runName}/comparison_shade/waterheight-$city-$subloc-$name_cont-$name_refl-shade\"\' \'container=\"$name_cont\"\' \'refl=\"$name_refl\"\' \'clouds=\"Observed\"\' \'site=\"$city-$subloc\"\' plots/plot_comparison_shade.ncl";
				}
			}
		}
	}

	if ($plot)
	{
		# Plotting - Container experiments
		print "## Running container experiment plots...\n";
		for (my $j = 0; $j < $exp_refl_size; $j++)
		{
			my $name_refl = $names_refl[$j];
			for (my $k = 0; $k < $exp_shade_size; $k++)
			{
				my $name_shade = $names_shade[$k];
				if ($cflag eq 'F')
				{
					for (my $l = 0; $l < $exp_cloud_size; $l++)
					{
						my $name_cloud = $names_cloud[$l];
						system "ncl -Q \'infile1=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-Coffee-$name_refl-$name_shade-$name_cloud.txt\"\' \'infile2=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-Bucket-$name_refl-$name_shade-$name_cloud.txt\"\' \'infile3=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-Barrel-$name_refl-$name_shade-$name_cloud.txt\"\' \'outfile1=\"$path_out_plots/${runName}/comparison_containers/temperature-$city-$subloc-$name_refl-$name_shade-$name_cloud-containers\"\' \'outfile2=\"$path_out_plots/${runName}/comparison_containers/waterheight-$city-$subloc-$name_refl-$name_shade-$name_cloud-containers\"\' \'refl=\"$name_refl\"\' \'shade=\"$name_shade\"\' \'clouds=\"$name_cloud\"\' \'site=\"$city-$subloc\"\' plots/plot_comparison_container.ncl";
					}
				}
				else
				{
					system "ncl -Q \'infile1=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-Coffee-$name_refl-$name_shade.txt\"\' \'infile2=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-Bucket-$name_refl-$name_shade.txt\"\' \'infile3=\"$path_out_data/${dataSrc}$city-$subloc/container_output_${dataSrc}$city-$subloc-Barrel-$name_refl-$name_shade.txt\"\' \'outfile1=\"$path_out_plots/${runName}/comparison_containers/temperature-$city-$subloc-$name_refl-$name_shade-containers\"\' \'outfile2=\"$path_out_plots/${runName}/comparison_containers/waterheight-$city-$subloc-$name_refl-$name_shade-containers\"\' \'refl=\"$name_refl\"\' \'shade=\"$name_shade\"\' \'clouds=\"Observed\"\' \'site=\"$city-$subloc\"\' plots/plot_comparison_container.ncl";
				}
			}
		}

		if ($cflag eq 'T')
        {
            my $ncols = 12;
            if ($sflag eq 'F') { $ncols = 11 }
            # Plotting - Cloud fraction and precipitation
            print "## Running cloud fraction and precipitation plots...\n";
            system "ncl -Q \'infile=\"input/container_input_${dataSrc}$city-$subloc.txt\"\' \'outfile=\"$path_out_plots/${runName}/clouds_precip/clouds-precip-$city-$subloc\"\' \'site=\"$city-$subloc\"\' \'ncols=$ncols\'  plots/plot_clouds_precip.ncl";
            # Plotting - Temperature and RH
            print "## Running temperature and RH plots...\n";
            system "ncl -Q \'infile=\"input/container_input_${dataSrc}$city-$subloc.txt\"\' \'outfile=\"$path_out_plots/${runName}/temp_and_rh/temp-rh-$city-$subloc\"\' \'site=\"$city-$subloc\"\' \'ncols=$ncols\'  plots/plot_temp_rh.ncl";
        }
        else
        {
            my $ncols = 9;
            if ($sflag eq 'F') { $ncols = 8 }
            # Plotting - Temperature and RH
            print "## Running temperature and RH plots...\n";
            system "ncl -Q \'infile=\"input/container_input_${dataSrc}$city-$subloc.txt\"\' \'outfile=\"$path_out_plots/${runName}/temp_and_rh/temp-rh-$city-$subloc\"\' \'site=\"$city-$subloc\"\' \'ncols=$ncols\'  plots/plot_temp_rh.ncl";
        }
	}
}

if ($plot)
{
	# Plotting - City experiments
	print "# Running city experiment plots...\n";
	for (my $i = 0; $i < $exp_cont_size; $i++)
	{
		my $name_cont = $names_cont[$i];
		for (my $j = 0; $j < $exp_refl_size; $j++)
		{
			my $name_refl = $names_refl[$j];
			for (my $k = 0; $k < $exp_shade_size; $k++)
			{
				my $name_shade = $names_shade[$k];
				if ($cflag eq 'F')
				{
					for (my $l = 0; $l < $exp_cloud_size; $l++)
					{
						my $name_cloud = $names_cloud[$l];
						system "ncl -Q \'infile1=\"$path_out_data/$cities[0]-$sublocs[0]/container_output_${dataSrc}$cities[0]-$sublocs[0]-$name_cont-$name_refl-$name_shade-$name_cloud.txt\"\' \'infile2=\"$path_out_data/$cities[1]-$sublocs[1]/container_output_${dataSrc}$cities[1]-$sublocs[1]-$name_cont-$name_refl-$name_shade-$name_cloud.txt\"\' \'infile3=\"$path_out_data/$cities[2]-$sublocs[2]/container_output_${dataSrc}$cities[2]-$sublocs[2]-$name_cont-$name_refl-$name_shade-$name_cloud.txt\"\' \'outfile1=\"$path_out_plots/${runName}/comparison_cities/temperature-$name_cont-$name_refl-$name_shade-$name_cloud\"\' \'outfile2=\"$path_out_plots/${runName}/comparison_cities/waterheight-$name_cont-$name_refl-$name_shade-$name_cloud\"\' \'shade=\"$name_shade\"\' \'refl=\"$name_refl\"\' \'clouds=\"$name_cloud\"\' \'container=\"$name_cont\"\' plots/plot_comparison_city.ncl";
					}
				}
				else
				{
					system "ncl -Q \'infile1=\"$path_out_data/$cities[0]-$sublocs[0]/container_output_${dataSrc}$cities[0]-$sublocs[0]-$name_cont-$name_refl-$name_shade.txt\"\' \'infile2=\"$path_out_data/$cities[1]-$sublocs[1]/container_output_${dataSrc}$cities[1]-$sublocs[1]-$name_cont-$name_refl-$name_shade.txt\"\' \'infile3=\"$path_out_data/$cities[2]-$sublocs[2]/container_output_${dataSrc}$cities[2]-$sublocs[2]-$name_cont-$name_refl-$name_shade.txt\"\' \'outfile1=\"$path_out_plots/${runName}/comparison_cities/temperature-$name_cont-$name_refl-$name_shade\"\' \'outfile2=\"$path_out_plots/${runName}/comparison_cities/waterheight-$name_cont-$name_refl-$name_shade\"\' \'shade=\"$name_shade\"\' \'refl=\"$name_refl\"\' \'clouds=\"Observed\"\' \'container=\"$name_cont\"\' plots/plot_comparison_city.ncl";
				}
			}
		}
	}
}
close MTOUTFILE;
close MWOUTFILE;

# Now open the summary statistics files to run some additional plots
print "# Writing summary statistics files and plots...\n";
my %data;
my @tvals = qw/ TMEAN TMAX TMIN CMEAN CMAX CMIN /;
open MTINFILE, "<stats/stats_temperature.txt" or die "Cannot open input stats temperature file: $!";
open MWINFILE, "<stats/stats_waterheight.txt" or die "Cannot open input stats water height file: $!";
my $junk = <MTINFILE>;
$junk = <MWINFILE>;
while (my $tline = <MTINFILE>)
{
	my $wline = <MWINFILE>;
	my @tparts = split /\s+/, $tline;
	my @wparts = split /\s+/, $wline;
	if ($cflag eq 'F')
	{
		$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{$tparts[4]}->{'TMEAN'} = $tparts[5];
		$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{$tparts[4]}->{'TMAX'} = $tparts[6];
		$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{$tparts[4]}->{'TMIN'} = $tparts[7];
		$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{$tparts[4]}->{'CMEAN'} = $tparts[8];
        $data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{$tparts[4]}->{'CMAX'} = $tparts[9];
        $data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{$tparts[4]}->{'CMIN'} = $tparts[10];
		$data{$wparts[0]}->{$wparts[1]}->{$wparts[2]}->{$wparts[3]}->{$tparts[4]}->{'WH'} = $wparts[5];
	}
	else
	{
		$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{'TMEAN'} = $tparts[4];
		$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{'TMAX'} = $tparts[5];
		$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{'TMIN'} = $tparts[6];
		$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{'CMEAN'} = $tparts[7];
        $data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{'CMAX'} = $tparts[8];
        $data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{'CMIN'} = $tparts[9];
		$data{$wparts[0]}->{$wparts[1]}->{$wparts[2]}->{$wparts[3]}->{'WH'} = $wparts[4];
	}
}
close MTINFILE;
close MWINFILE;

foreach my $name_refl (@names_refl)
{
	foreach my $name_shade (@names_shade)
	{
		open WTEMPFILE, ">temp_plot_wh.txt" or die "Cannot open temp plot file: $!";
		foreach my $tval (@tvals)
		{
			open TTEMPFILE, ">temp_plot_t_$tval.txt" or die "Cannot open temp plot file: $!";
			foreach my $name_cont (@names_cont)
			{
				if ($cflag eq 'F')
				{
					foreach my $name_cloud (@names_cloud)
					{
						printf TTEMPFILE "%6.2f %6.2f %6.2f ", $data{"$cities[0]-$sublocs[0]"}->{$name_cont}->{$name_refl}->{$name_shade}->{$name_cloud}->{$tval}, $data{"$cities[1]-$sublocs[1]"}->{$name_cont}->{$name_refl}->{$name_shade}->{$name_cloud}->{$tval}, $data{"$cities[2]-$sublocs[2]"}->{$name_cont}->{$name_refl}->{$name_shade}->{$name_cloud}->{$tval};
					}
				}
				else
				{
					printf TTEMPFILE "%6.2f %6.2f %6.2f ", $data{"$cities[0]-$sublocs[0]"}->{$name_cont}->{$name_refl}->{$name_shade}->{$tval}, $data{"$cities[1]-$sublocs[1]"}->{$name_cont}->{$name_refl}->{$name_shade}->{$tval}, $data{"$cities[2]-$sublocs[2]"}->{$name_cont}->{$name_refl}->{$name_shade}->{$tval};
				}
			}
			print TTEMPFILE "\n";
			close TTEMPFILE;

			# Plotting - Temperature bar plots
			if ($plot)
			{
				if ($cflag eq 'F')
				{
					foreach my $name_cloud (@names_cloud)
					{
						system "ncl -Q \'infile=\"temp_plot_t_$tval.txt\"\' \'outfile=\"$path_out_plots/${runName}/bar_temperature/$tval-$name_refl-$name_shade-$name_cloud\"\' \'shade=\"$name_shade\"\' \'refl=\"$name_refl\"\' \'clouds=\"$name_cloud\"\' plots/plot_bar_t.ncl";
					}
				}
				else
				{
					system "ncl -Q \'infile=\"temp_plot_t_$tval.txt\"\' \'outfile=\"$path_out_plots/${runName}/bar_temperature/$tval-$name_refl-$name_shade\"\' \'shade=\"$name_shade\"\' \'refl=\"$name_refl\"\' \'clouds=\"Observed\"\' plots/plot_bar_t.ncl";
				}
			}
		}

		# Plotting - Summary DTR
		if ($plot)
		{
			if ($cflag eq 'F')
			{
				foreach my $name_cloud (@names_cloud)
				{
					system "ncl -Q \'infile1=\"temp_plot_t_TMAX.txt\"\' \'infile2=\"temp_plot_t_TMIN.txt\"\' \'outfile=\"$path_out_plots/${runName}/bar_dtr/DTR-$name_refl-$name_shade-$name_cloud\"\' \'shade=\"$name_shade\"\' \'refl=\"$name_refl\"\' \'clouds=\"$name_cloud\"\' plots/plot_bar_dtr.ncl";
				}
			}
			else
			{
				system "ncl -Q \'infile1=\"temp_plot_t_TMAX.txt\"\' \'infile2=\"temp_plot_t_TMIN.txt\"\' \'outfile=\"$path_out_plots/${runName}/bar_dtr/DTR-$name_refl-$name_shade\"\' \'shade=\"$name_shade\"\' \'refl=\"$name_refl\"\' \'clouds=\"Observed\"\' plots/plot_bar_dtr.ncl";
			}
		}

		foreach my $name_cont (@names_cont)
		{
			if ($cflag eq 'F')
			{
				foreach my $name_cloud (@names_cloud)
				{
					printf WTEMPFILE "%6.2f %6.2f %6.2f ", $data{"$cities[0]-$sublocs[0]"}->{$name_cont}->{$name_refl}->{$name_shade}->{$name_cloud}->{'WH'}, $data{"$cities[1]-$sublocs[1]"}->{$name_cont}->{$name_refl}->{$name_shade}->{$name_cloud}->{'WH'}, $data{"$cities[2]-$sublocs[2]"}->{$name_cont}->{$name_refl}->{$name_shade}->{$name_cloud}->{'WH'};
				}
			}
			else
			{
				printf WTEMPFILE "%6.2f %6.2f %6.2f ", $data{"$cities[0]-$sublocs[0]"}->{$name_cont}->{$name_refl}->{$name_shade}->{'WH'}, $data{"$cities[1]-$sublocs[1]"}->{$name_cont}->{$name_refl}->{$name_shade}->{'WH'}, $data{"$cities[2]-$sublocs[2]"}->{$name_cont}->{$name_refl}->{$name_shade}->{'WH'};
			}

			# Plotting - Summary frequency distribution of temperature by city
			if ($plot)
			{
				if ($cflag eq 'F')
				{
					foreach my $name_cloud (@names_cloud)
					{
						system "ncl -Q \'infile1=\"$path_out_data/$cities[0]-$sublocs[0]/stats_$cities[0]-$sublocs[0]-$name_cont-$name_refl-$name_shade-$name_cloud.txt\"\' \'infile2=\"$path_out_data/$cities[1]-$sublocs[1]/stats_$cities[1]-$sublocs[1]-$name_cont-$name_refl-$name_shade-$name_cloud.txt\"\' \'infile3=\"$path_out_data/$cities[2]-$sublocs[2]/stats_$cities[2]-$sublocs[2]-$name_cont-$name_refl-$name_shade-$name_cloud.txt\"\' \'outfile=\"$path_out_plots/${runName}/freq_distribution/tmean-$name_cont-$name_refl-$name_shade-$name_cloud\"\' \'shade=\"$name_shade\"\' \'refl=\"$name_refl\"\' \'clouds=\"$name_cloud\"\' \'container=\"$name_cont\"\' plots/plot_freq_dist_t.ncl";
					}
				}
				else
				{
					system "ncl -Q \'infile1=\"$path_out_data/$cities[0]-$sublocs[0]/stats_$cities[0]-$sublocs[0]-$name_cont-$name_refl-$name_shade.txt\"\' \'infile2=\"$path_out_data/$cities[1]-$sublocs[1]/stats_$cities[1]-$sublocs[1]-$name_cont-$name_refl-$name_shade.txt\"\' \'infile3=\"$path_out_data/$cities[2]-$sublocs[2]/stats_$cities[2]-$sublocs[2]-$name_cont-$name_refl-$name_shade.txt\"\' \'outfile=\"$path_out_plots/${runName}/freq_distribution/tmean-$name_cont-$name_refl-$name_shade\"\' \'shade=\"$name_shade\"\' \'refl=\"$name_refl\"\' \'clouds=\"Observed\"\' \'container=\"$name_cont\"\' plots/plot_freq_dist_t.ncl";
				}
			}
		}
		print WTEMPFILE "\n";
		close WTEMPFILE;
		# Plotting - Summary water height distribution
		if ($plot)
		{
			if ($cflag eq 'F')
			{
				foreach my $name_cloud (@names_cloud)
				{
					system "ncl -Q \'infile=\"temp_plot_wh.txt\"\' \'outfile=\"$path_out_plots/${runName}/bar_waterheight/wh-$name_refl-$name_shade-$name_cloud\"\' \'shade=\"$name_shade\"\' \'refl=\"$name_refl\"\' \'clouds=\"$name_cloud\"\' plots/plot_bar_wh.ncl";
				}
			}
			else
			{
				system "ncl -Q \'infile=\"temp_plot_wh.txt\"\' \'outfile=\"$path_out_plots/${runName}/bar_waterheight/wh-$name_refl-$name_shade\"\' \'shade=\"$name_shade\"\' \'refl=\"$name_refl\"\' \'clouds=\"Observed\"\' plots/plot_bar_wh.ncl";
			}
		}
	}
}

# Clean up
print "# Cleaning up...\n";
system "rm temp_plot_t_*.txt temp_plot_wh.txt param_file.txt";
print "# Finished!\n";

# Subroutines
sub daily_energy_balance_stats
{
	my ($einfile, $eoutfile, $soutfile, $citysubloc, $name_cont, $name_refl, $name_shade, $cflag, $name_cloud) = @_;
	# Construct average daily energy balance and some other statistics
	my %ebal;
	my %tw;
	my %tc;
	my %ta;
	my %evap;
	my $whc = 0;
	open EINFILE, "<$einfile" or die "Cannot open input file for energy balance: $!";
	for (my $x = 0; $x < 7; $x++)
	{
		my $junk = <EINFILE>; # Header lines
	}
	while (my $line = <EINFILE>)
	{
		chomp ($line);
		my @parts = split /\s+/, $line;
		my $md = $parts[1].$parts[2];
		my $hour = $parts[3];
		if ($parts[24] != $miss)
		{
			push (@{$ebal{$hour}->{'SWDOWN'}}, $parts[4]);
			push (@{$ebal{$hour}->{'SWSIDE'}}, $parts[5]);
			push (@{$ebal{$hour}->{'LWDOWN'}}, $parts[6]);
			push (@{$ebal{$hour}->{'LWUP'}}, $parts[7]);
			push (@{$ebal{$hour}->{'LWIN'}}, $parts[8]);
			push (@{$ebal{$hour}->{'LWOUT'}}, $parts[9]);
			push (@{$ebal{$hour}->{'SHUP'}}, $parts[10]);
			push (@{$ebal{$hour}->{'SHOUT'}}, $parts[11]);
			push (@{$ebal{$hour}->{'LHUP'}}, $parts[12]);
			push (@{$ebal{$hour}->{'GHW'}}, $parts[13]);
			push (@{$ebal{$hour}->{'GHC'}}, $parts[14]);
			push (@{$ebal{$hour}->{'COND'}}, $parts[15]);
			push (@{$ebal{$hour}->{'BAL_W'}}, $parts[16]);
			push (@{$ebal{$hour}->{'BAL_C'}}, $parts[17]);
			my $sw = $parts[4] + $parts[5];
			push (@{$ebal{$hour}->{'SW'}}, $sw);
			my $lw = $parts[6] - $parts[7] + $parts[8] - $parts[9];
			push (@{$ebal{$hour}->{'LW'}}, $lw);
			my $sh = -1.*($parts[10] + $parts[11]);
			push (@{$ebal{$hour}->{'SH'}}, $sh);
			push (@{$evap{$md}}, $parts[19]);
			push (@{$ebal{$hour}->{'TG'}}, $parts[21]);
			push (@{$ebal{$hour}->{'TA'}}, $parts[22]);
			push (@{$ebal{$hour}->{'TCON'}}, $parts[23]);
			push (@{$ebal{$hour}->{'TW'}}, $parts[24]);
			push (@{$tw{$md}}, $parts[24]);
			push (@{$tc{$md}}, $parts[23]);
			push (@{$ta{$md}}, $parts[22]);
			if ($parts[20] >= 0.015)
			{
				$whc++;
			}
		}
	}
	close EINFILE;

	# Write out daily average energy balance file
	open EOUTFILE, ">$eoutfile" or die "Cannot open output file for energy balance: $!";
	print EOUTFILE "Hour,SWDOWN,SWSIDE,LWDOWN,LWUP,LWIN,LWOUT,SHUP,SHOUT,LHUP,GHW,GHC,COND,BAL_W,BAL_C,SW,LW,SH,TG,TA,TCON,TW\n";
	for (my $hr = 0; $hr < 24; $hr++)
	{
		printf EOUTFILE "%02s ", $hr;
		my @terms = qw/ SWDOWN SWSIDE LWDOWN LWUP LWIN LWOUT SHUP SHOUT LHUP GHW GHC COND BAL_W BAL_C SW LW SH TG TA TCON TW /;
		foreach my $term (@terms)
		{
			my $estats = Statistics::Descriptive::Full->new();
			$estats->add_data(@{$ebal{$hr}->{$term}});
			my $mean = $estats->mean();
			printf EOUTFILE "%8.2f ", $mean;
		}
		print EOUTFILE "\n";
	}
	close EOUTFILE;

	# Calculate some basic statistics and output to file
	open SOUTFILE, ">$soutfile" or die "Cannot open output file for stats: $!";
	my %tastats;
	$tastats{'TMEAN'} = Statistics::Descriptive::Full->new();
	$tastats{'TMAX'} = Statistics::Descriptive::Full->new();
	$tastats{'TMIN'} = Statistics::Descriptive::Full->new();
	$tastats{'CMEAN'} = Statistics::Descriptive::Full->new();
    $tastats{'CMAX'} = Statistics::Descriptive::Full->new();
    $tastats{'CMIN'} = Statistics::Descriptive::Full->new();
	$tastats{'TAMEAN'} = Statistics::Descriptive::Full->new();
	$tastats{'EVAP'} = Statistics::Descriptive::Full->new();

	foreach my $md (sort keys %tw)
	{
		my $tstats = Statistics::Descriptive::Full->new();
		my $astats = Statistics::Descriptive::Full->new();
		my $cstats = Statistics::Descriptive::Full->new();
		$tstats->add_data(@{$tw{$md}});
		$astats->add_data(@{$ta{$md}});
		$cstats->add_data(@{$tc{$md}});
		$tastats{'EVAP'}->add_data(@{$evap{$md}});
		my $tmean = $tstats->mean();
		my $tmax = $tstats->max();
		my $tmin = $tstats->min();
		my $tamean = $astats->mean();
		my $cmean = $cstats->mean();
        my $cmax = $cstats->max();
        my $cmin = $cstats->min();
		if (scalar(@{$tw{$md}}) == 24)
		{
			$tastats{'TMEAN'}->add_data($tmean);
			$tastats{'TMAX'}->add_data($tmax);
			$tastats{'TMIN'}->add_data($tmin);
			$tastats{'CMEAN'}->add_data($cmean);
            $tastats{'CMAX'}->add_data($cmax);
            $tastats{'CMIN'}->add_data($cmin);
			$tastats{'TAMEAN'}->add_data($tamean);
		}
	}
	print SOUTFILE "Average Water Daily Temperature (Mean, Max, Min): \n";
	my $tmean = $tastats{'TMEAN'}->mean();
	my $tmax = $tastats{'TMAX'}->mean();
	my $tmin = $tastats{'TMIN'}->mean();
	my $cmean = $tastats{'CMEAN'}->mean();
    my $cmax = $tastats{'CMAX'}->mean();
    my $cmin = $tastats{'CMIN'}->mean();
	my $tamean = $tastats{'TAMEAN'}->mean();
	my $evapmean = $tastats{'EVAP'}->mean();
	printf SOUTFILE "%6.2f %6.2f %6.2f\n", $tmean, $tmax, $tmin;
	print SOUTFILE "Average Container Daily Temperature (Mean, Max, Min): \n";
	printf SOUTFILE "%6.2f %6.2f %6.2f\n", $cmean, $cmax, $cmin;
	if ($cflag eq 'F')
	{
		printf MTOUTFILE "%-30s %-10s %-15s %-10s %-15s %-6.2f %-6.2f %-6.2f %-6.2f%-6.2f %-6.2f %-6.2f %-6.4f\n", $citysubloc, $name_cont, $name_refl, $name_shade, $name_cloud, $tmean, $tmax, $tmin, $cmean, $cmax, $cmin, $tamean, $evapmean;
	}
	else
	{
		printf MTOUTFILE "%-30s %-10s %-15s %-10s %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.4f\n", $citysubloc, $name_cont, $name_refl, $name_shade, $tmean, $tmax, $tmin, $cmean, $cmax, $cmin, $tamean, $evapmean;
	}
	print SOUTFILE "Number of days water level > 15 mm: \n";
	$whc = $whc / 24;
	printf SOUTFILE "%6.2f\n", $whc;
	if ($cflag eq 'F')
	{
		printf MWOUTFILE "%-30s %-10s %-15s %-10s %-15s %-6.2f\n", $citysubloc, $name_cont, $name_refl, $name_shade, $name_cloud, $whc;
	}
	else
	{
		printf MWOUTFILE "%-30s %-10s %-15s %-10s %-6.2f\n", $citysubloc, $name_cont, $name_refl, $name_shade, $whc;
	}
	print SOUTFILE "Frequency Distribution: \n";
	my @bins = (10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50); # new bins
	my $fd = $tastats{'TMEAN'}->frequency_distribution_ref(\@bins);
	for (sort {$a <=> $b} keys %{$fd})
	{
		print SOUTFILE "$_  $fd->{$_}\n";
	}
	close SOUTFILE;
}
