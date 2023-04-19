Tutorial
=========

.. code-block:: python

	from cfdna import cfDNA, cfDNA_summary
	from cfdna.Plot import plot_plt
	import numpy as np
	
	# Read sam file
	f = cfDNA("test_bam.bam", verbose=True)
	
	# Determine CNVs
	cnvs = f.call_cnvs()
	
	# Determine window protection score (WPS)
	wps = f.wps()
	
	# Summarize
	summary = f.summarize()
	
	# Write summary
	summary.write_h5("test.h5")
	
	# Read summary
	summary = cfDNA_summary("test.h5")
	
	# Plot
	plot_plt.plot_summary(summary)