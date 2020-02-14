function springs_setup()

dir_solver = fullfile('.','springs_solver') ;
dir_solver_obj = fullfile(dir_solver,'obj') ;
if ~exist(dir_solver_obj)
	mkdir(dir_solver_obj)
end
system(sprintf('make -C %s clean',dir_solver)) ;
system(sprintf('make -C %s',dir_solver)) ;

end