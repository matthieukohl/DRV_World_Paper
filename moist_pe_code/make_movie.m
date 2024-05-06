% make video

video = VideoWriter('movie_pe_r001_eps04_final_paper.avi');
video.FrameRate = 10;
open(video)

test = dir([pwd,'/images_movie']);

for ii = 1:320
    I = imread([pwd,'/images_movie/frame_',num2str(ii),'.png']);
    writeVideo(video,I);
end

close(video)
    