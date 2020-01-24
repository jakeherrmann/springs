function [ nodes , springs ] = springs_solve( nodes , springs )

dir_input  = 'springs_input' ;
dir_output = 'springs_output' ;
precision = 'double' ;

springs_write( dir_input , precision , nodes , springs ) ;
exe_springs_solver = fullfile('.','springs_solver','springs_solver.exe') ;
sys_command = sprintf( '%s %s %s' , exe_springs_solver , dir_input , dir_output ) ;
system(sys_command) ;
[ nodes , springs ] = springs_read( dir_output ) ;

end