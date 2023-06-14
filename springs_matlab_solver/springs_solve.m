function [ nodes , springs ] = springs_solve( nodes , springs , varargin )

% default/specified solver parameters
options = springs_default_options() ;
if numel(varargin) > 0
	for ff = fieldnames(varargin{1})'
		options.(ff{1}) = varargin{1}.(ff{1}) ;
	end
end
if numel(varargin) > 1
	use_parallel = varargin{2} ;
end
if numel(varargin) > 2
	verbose = varargin{3} ;
end

exe_springs_solver = fullfile('.','springs_solver','springs_solver.exe') ;
if ~exist( exe_springs_solver ,'file')
	fprintf( 'ERROR: Compiled executable not found.  Attempting to compile.\n' )
	springs_setup()
end

dir_input  = 'springs_input' ;
dir_output = 'springs_output' ;

fids = fopen('all') ;
filenames = arrayfun( @fopen , fids ,'UniformOutput',false) ;
fids = fids( contains(filenames,{dir_input,dir_output}) ) ;
if ~isempty( fids )
	fclose( fids ) ;
end
if ~exist( dir_input ,'dir')
	mkdir( dir_input )
else
	delete(fullfile(dir_input,'*'))
end
if ~exist( dir_output ,'dir')
	mkdir( dir_output )
else
	delete(fullfile(dir_output,'*')) ;
end

springs_write( dir_input , nodes , springs , options ) ;
setenv( 'OMP_STACKSIZE'   , '64M'    ) ;
setenv( 'OMP_PLACES'      , 'cores'  ) ;
setenv( 'OMP_WAIT_POLICY' , 'active' ) ;
setenv( 'OMP_PROC_BIND'   , 'true'   ) ;
sys_command = sprintf( '%s %s %s --verbose %d --parallel %d' , exe_springs_solver , dir_input , dir_output , verbose , use_parallel ) ;
if verbose
	system(sys_command) ;
else
	[~,~] = system(sys_command) ;
end
[ nodes , springs ] = springs_read( dir_output ) ;

end