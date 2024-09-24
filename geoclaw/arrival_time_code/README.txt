GeoClaw Arrival Time Code:

setrun.py # example set run to run geoclaw tsunami model- specified wave arrival tolerance thresholds set in line 638.
		## fgmax_tools.py need to be modified with arrival_tols_h variable names to run successfully. 

process_fgmax_arrival_time_example.ipynb # jupyter notebook to convert results to netCDF format for plotting in arcGIS products.

plot_fgmax.py # post-processing script to generate example_plots in _output.

make_summary_Oct31.py # script to generate _output/gauges_report_timesarrival.txt summary report
