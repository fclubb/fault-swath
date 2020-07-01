#---------------------------------------------------------------------#
# Run plotting along San Andreas Fault
# Developed by Fiona Clubb
# Durham University
#---------------------------------------------------------------------#

# setting backend to run on server
import matplotlib
matplotlib.use('Agg')
import os
import sys
import pandas as pd
import swath_profile as swath

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to plot the analysis along the fault for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -dir flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("You also need the -fname flag which will give the prefix of the raster files.")
    print("For help type:")
    print("   python run-swath-profile.py -h\n")
    print("=======================================================================\n\n ")

if __name__ == '__main__':

    # If there are no arguments, send to the welcome screen
    if (len(sys.argv) < 1):
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()

    # The location of the data files
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error.")

    # options for analysis
    parser.add_argument("-so", "--stream_order", type=int, default=3, help="The stream order of your channel and hillslope profiles. Default = 3")

    # Different options for plotting
    parser.add_argument("-ch", "--channels", type=bool, default=False, help="Plot channel slope vs. distance along the fault")
    parser.add_argument("-b", "--basins", type=bool, default=False, help="Plot basin orientation along the fault.")
    parser.add_argument("-tc", "--thermochron", type=bool, default=False, help="If this is true I'll make plots of the channel slopes against uplift rates from thermochron")
    parser.add_argument("-gps", "--gps", type=bool, default=False, help="If this is true I'll make plots of the channel slopes vs. gps data")
    parser.add_argument("-insar", "--insar", type=bool, default=False, help="If this is true I'll make plots of the channel slopes vs. InSAR data")
    parser.add_argument("-lith", "--lithology", type=bool, default=False, help="If this is true I'll make plots of the channel slopes separated by lithology")
    parser.add_argument("-hs", "--hillslopes", type=bool, default=False, help="If this is true I'll make plots of hillslope gradient vs. distance")
    parser.add_argument("-multiple_so", "--multiple_so", type=bool, default=False, help="If this is true I'll plot all the stream orders")
    parser.add_argument("-eq", "--earthquakes", type=bool, default=False, help="If this is true I'll plot the distribution of earthquakes along the fault")
    parser.add_argument("-az", "--azimuth", type=bool, default=False, help="If this is true I'll plot the fault aziumth along strike.")

    args = parser.parse_args()

    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()
    else:
        fname_prefix = args.fname_prefix

    # get the base directory
    if args.base_directory:
        DataDirectory = args.base_directory
        # check if you remembered a / at the end of your path_name
        # if not DataDirectory.endswith("/"):
        #     print("You forgot the '/' at the end of the directory, appending...")
        #     DataDirectory = DataDirectory+"/"
    else:
        print("WARNING! You haven't supplied the data directory. I'm using the current working directory.")
        DataDirectory = os.getcwd()
    base_dir = '/home/bjdd72/TopographicData/TopographicData/san_andreas/'

    # print the arguments that you used to an output file for reproducibility
    with open(DataDirectory+args.fname_prefix+'_report.csv', 'w') as output:
        for arg in vars(args):
            output.write(str(arg)+','+str(getattr(args, arg))+'\n')
        output.close()

    # read in the profiles and check if you've calculated the distance along the fault
    profile_csv = DataDirectory+fname_prefix+'_profiles_SO{}.csv'.format(args.stream_order)
    output_csv=DataDirectory+fname_prefix+'_profiles_fault_dist_SO{}.csv'.format(args.stream_order)
    baseline_shapefile='SanAndreasFault.shp'
    output_shapefile='SanAndreasPoints.shp'

    # check if the fault dist csv already exists
    if not os.path.isfile(output_csv):
        points, distances = swath.get_points_along_line(DataDirectory,baseline_shapefile,output_shapefile,n=512)
        coeffs = swath.get_orthogonal_coefficients(points)
        profile_df = pd.read_csv(profile_csv)
        swath.bisection_method(points, coeffs, distances, profile_df, output_csv)

    # labels
    labels_csv=base_dir+'Uplift_rates/placenames.csv'
    labels_df = pd.read_csv(labels_csv)
    swath.get_distance_along_fault_from_points(DataDirectory, baseline_shapefile, labels_df, labels_csv)

    # horizontal slip rate
    slip_rate_csv = base_dir+'Slip_rates/Tong_2013_InSAR.csv'
    output_sr_csv = base_dir+'Slip_rates/Tong_2013_InSAR_fault_dist.csv'
    if not os.path.isfile(output_sr_csv):
        sr_df = pd.read_csv(slip_rate_csv)
        swath.get_distance_along_fault_from_points(DataDirectory, baseline_shapefile, sr_df, output_sr_csv)

    # seismic data
    eq_csv = base_dir+'Earthquakes/SAF_earthquakes_50km.csv'
    output_eq_csv =  base_dir+'Earthquakes/SAF_earthquakes_50km_dist.csv'
    if not os.path.isfile(output_eq_csv):
        eq_dr = pd.read_csv(eq_csv)
        swath.get_distance_along_fault_from_points(DataDirectory, baseline_shapefile, eq_df, output_eq_csv)

    # gps data
    gps_csv=base_dir+'Uplift_rates/gps/MIDAS_IGS08_SAF_50km.csv'
    output_gps_csv=base_dir+'Uplift_rates/gps/MIDAS_IGS08_SAF_50km_dist.csv'
    if not os.path.isfile(output_gps_csv):
        gps_df = pd.read_csv(gps_csv)
        swath.get_distance_along_fault_from_points(DataDirectory, baseline_shapefile, gps_df, output_gps_csv)


    # channel slope plotting
    if args.channels:
        swath.plot_channel_slopes_along_fault(DataDirectory, fname_prefix, args.stream_order, output_csv, labels_csv, output_sr_csv, output_eq_csv)

    if args.basins:
        basins = DataDirectory+fname_prefix+'_basins_WGS84.shp'
        swath.plot_basin_orientation_along_fault(DataDirectory, fname_prefix, basins, baseline_shapefile, output_shapefile, output_sr_csv)

    # hillslope plotting
    if args.hillslopes:
        hillslope_csv= DataDirectory+fname_prefix+'_hillslopes_SO{}.csv'.format(args.stream_order)
        output_hillslope_csv=DataDirectory+fname_prefix+'_hillslopes_SO{}_dist.csv'.format(args.stream_order)
        if not os.path.isfile(output_hillslope_csv):
            points, distances = swath.get_points_along_line(n=512)
            coeffs = swath.get_orthogonal_coefficients(points)
            hillslope_df = pd.read_csv(hillslope_csv)
            swath.bisection_method(points, coeffs, distances, hillslope_df, output_hillslope_csv)
        # do the plotting
        swath.plot_hillslopes_along_fault(output_hillslope_csv)

    # lithology
    if args.lithology:
        lithology_raster=base_dir+'Lithology/ca_geol_simple_utm.tif'
        lithology_shp=base_dir+'Lithology/ca_geol_simple_dissolved.shp'
        output_lith_csv = DataDirectory+fname_prefix+'_profiles_lithology_SO{}.csv'.format(args.stream_order)
        if not os.path.isfile(output_lith_csv):
            swath.burn_lithology_to_river_df(output_csv, output_lith_csv, lithology_raster)
        # plotting
        #swath.plot_lithology_shapefile(DataDirectory,lithology_shp,baseline_shapefile)
        swath.plot_channel_slopes_uniform_lithology(DataDirectory, fname_prefix, output_lith_csv, labels_csv, args.stream_order)
        #swath.plot_lithology_deltas(DataDirectory, fname_prefix, output_lith_csv, labels_csv, args.stream_order)

    # multiple stream orders
    if args.multiple_so:
        swath.plot_channel_slopes_multiple_SO(DataDirectory,fname_prefix,labels_csv)

    # seismic data
    #if args.eq:
    #    swath.plot_earthquakes_along_fault(DataDirectory, fname_prefix, output_eq_csv)


    print("Done, enjoy your plots!")
