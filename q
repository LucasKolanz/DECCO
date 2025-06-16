[1mdiff --git a/data/DECCOData.cpp b/data/DECCOData.cpp[m
[1mindex f64377c..77d1d15 100644[m
[1m--- a/data/DECCOData.cpp[m
[1m+++ b/data/DECCOData.cpp[m
[36m@@ -1352,6 +1352,7 @@[m [mbool DECCOData::write_checkpoint()[m
 		checkpt_file = filename.substr(0,filename.find_last_of('_')+1)+"checkpoint.txt";[m
 	}[m
 [m
[32m+[m	[32m// std::cerr<<"CREATING FILE: "<<checkpt_file<<std::endl;[m
 	std::ofstream output(checkpt_file);[m
 	if (not output.good()) [m
 	{[m
