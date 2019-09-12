
dir_image = '/Users/jake/Documents/GitHub/springNet/springs/TEST_RESULT/' ;
file_list = dir(fullfile(dir_image,'ITER*.png')) ;
clear img
for ii = numel(file_list) : -1 : 1
	img{ii} = imread( fullfile(dir_image,file_list(ii).name) ) ;
end

%%

vid = VideoWriter( fullfile(dir_image,'test.mp4') , 'MPEG-4' ) ;
vid.FrameRate = 30 ;
vid.open()
for ii = 1 : numel(img)
	rep = 1 ;
% 	if     ii>255; rep = 1 ;
% 	elseif ii>250; rep = 2 ;
% 	elseif ii>245; rep = 3 ;
% 	elseif ii>235; rep = 4 ;
% 	elseif ii>230; rep = 3 ;
% 	elseif ii>225; rep = 2 ;
% 	else           rep = 1 ;
% 	end
	for nn = 1 : rep
		vid.writeVideo( img{ii} )
	end
end
vid.close()