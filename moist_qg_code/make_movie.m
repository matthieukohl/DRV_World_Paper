% make video

video = VideoWriter('movie_qg_r001_linear_tau_damping_high_temp_regular_output.avi');
video.FrameRate = 8;
open(video)

test = dir([pwd,'/images_movie']);

for ii = 1:240
    I = imread(['images_movie/frame',num2str(ii),'.png']);
    writeVideo(video,I);
end

close(video)
    