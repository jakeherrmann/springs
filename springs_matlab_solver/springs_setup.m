function springs_setup()

if ~exist(fullfile('springs_solver','obj'),'dir')
	mkdir(fullfile('springs_solver','obj'))
end
system('make -C ./springs_solver clean') ;
system('make -C ./springs_solver') ;

end