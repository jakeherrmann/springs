function [ nodes , springs ] = springs_solve( nodes , springs , varargin )

% default/specified solver parameters
options = springs_default_options() ;
if numel(varargin) > 0
	for ff = fieldnames(varargin{1})'
		options.(ff{1}) = varargin{1}.(ff{1}) ;
	end
end
dir_input  = 'springs_input' ;
dir_output = 'springs_output' ;
if numel(varargin) > 1
	dir_input = varargin{2} ;
end
if numel(varargin) > 2
	dir_output = varargin{3} ;
end

exe_springs_solver = fullfile('.','springs_solver','springs_solver.exe') ;
if ~exist( exe_springs_solver ,'file')
	fprintf( 'ERROR: Compiled executable not found.  Attempting to compile.\n' )
	springs_setup()
end

check_files = {
	'network_parameters.txt'
	'network_nodes.dat'
	'network_springs.dat'
	} ;
fids = fopen('all') ;
filenames = arrayfun( @fopen , fids ,'UniformOutput',false) ;
fids = fids( contains(filenames,{dir_input,dir_output}) ) ;
if ~isempty( fids )
	fclose( fids ) ;
end
if ~exist( dir_input ,'dir')
	mkdir( dir_input )
else
	rm_files = fullfile(dir_input,check_files) ;
	delete( rm_files{:} )
end
if ~exist( dir_output ,'dir')
	mkdir( dir_output )
else
	rm_files = fullfile(dir_output,check_files) ;
	delete( rm_files{:} )
end

springs_write( dir_input , nodes , springs , options ) ;
setenv( 'OMP_STACKSIZE'   , '64M'    ) ;
setenv( 'OMP_PLACES'      , 'cores'  ) ;
setenv( 'OMP_WAIT_POLICY' , 'active' ) ;
setenv( 'OMP_PROC_BIND'   , 'true'   ) ;
sys_command = sprintf( '%s %s %s --verbose %d --parallel %d' , exe_springs_solver , dir_input , dir_output , options.verbose , options.use_parallel ) ;
if options.verbose
	sys_failure = system( sys_command ) ;
else
	[ sys_failure , ~ ] = system( sys_command ) ;
end
if sys_failure
	fprintf( 'SYSTEM() ERROR CODE: %d\n' , sys_failure )
	if ispc()
		fprintf( 'WIN64: attempting call from .bat file\n' )
		ml_path = fullfile( matlabroot() , 'bin' , 'win64' ) ;
		batfile = 'dev_solver.bat' ;
		fid = fopen(batfile,'w') ;
		fprintf(fid,'set path=%%path:%s;=%%\n',ml_path) ;
		fprintf(fid,'%s\n',sys_command) ;
		fclose( fid ) ;
		if options.verbose
			system( batfile ) ;
		else
			[ ~ , ~ ] = system( batfile ) ;
		end
	end
end

[ nodes , springs ] = springs_read( dir_output ) ;

end