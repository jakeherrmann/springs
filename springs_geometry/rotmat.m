function [ R ] = rotmat( rotvec )

t = rotvec(1) ;
Rx = [
	1 0 0
	0 cos(t) -sin(t)
	0 sin(t) cos(t)
	] ;

t = rotvec(2) ;
Ry = [
	cos(t) 0 sin(t)
	0 1 0
	-sin(t) 0 cos(t)
	] ;

t = rotvec(3) ;
Rz = [
	cos(t) -sin(t) 0
	sin(t) cos(t) 0
	0 0 1
	] ;

R = Rz * Ry * Rx ;

end