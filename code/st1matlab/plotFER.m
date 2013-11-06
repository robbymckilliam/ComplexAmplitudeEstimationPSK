
% load('results/st1_5p_da_20130909T114133.mat')
% figure
% semilogy(EbNodB_vec,numFrameErr./numFramesSimulated,'rx-','DisplayName','data-aided, 5 pilots')
% hold on, grid on
% load('results/st1_10p_da_20130909T114110.mat')
% semilogy(EbNodB_vec,numFrameErr./numFramesSimulated,'rx-','DisplayName','data-aided, 10 pilots')
% load('results/st1_15p_da_20130909T114048.mat')
% semilogy(EbNodB_vec,numFrameErr./numFramesSimulated,'rx-','DisplayName','data-aided, 15 pilots')
% load('results/st1_5p_Mackenthun_20130909T110106.mat')
% semilogy(EbNodB_vec,numFrameErr./numFramesSimulated,'go-','DisplayName','Mackenthun, 5 pilots')
% load('results/st1_10p_Mackenthun_20130909T110215.mat')
% semilogy(EbNodB_vec,numFrameErr./numFramesSimulated,'go-','DisplayName','Mackenthun, 10 pilots')
% load('results/st1_15p_Mackenthun_20130909T110746.mat')
% semilogy(EbNodB_vec,numFrameErr./numFramesSimulated,'go-','DisplayName','Mackenthun, 15 pilots')
% load('results/st1_5p_pc_20130909T151318.mat')
% semilogy(EbNodB_vec,numFrameErr./numFramesSimulated,'ks-','DisplayName','perfect channel')
% 
% xlabel('Eb/No [dB]')
% ylabel('FER')

%%
clear

% dirname = 'resultsOLD/'; fileCountExpected = 10;
% dirname = 'results/'; fileCountExpected = 10;
dirname = 'resultsCOMBINED/'; fileCountExpected = 20;

% basefilename = {'st1_5p_da_2013', 'st1_10p_da_2013', 'st1_15p_da_2013', ...
%   'st1_5p_Mackenthun_2013', 'st1_10p_Mackenthun_2013', 'st1_15p_Mackenthun_2013', ...
%   'st1_5p_pc_2013'};

basefilename = {'st1_5p_da_job', 'st1_10p_da_job', 'st1_15p_da_job', ...
  'st1_5p_Mackenthun_job', 'st1_10p_Mackenthun_job', 'st1_15p_Mackenthun_job', ...
  'st1_5p_pc_job'};
plotstyle    = {'r-x', 'r-x', 'r-x', ...
  'g-o', 'g-o', 'g-o', ...
  'k-d'};

list = what(dirname);
figure
for bsn_idx = 1:length(basefilename)
  numFrameErr_cum         = 0;
  numFramesSimulated_cum  = 0;
  fileCount = 0;
  for idx = 1:length(list.mat)
    k = strfind(list.mat{idx},basefilename{bsn_idx});
    if ~isempty(k)
      fileCount = fileCount + 1;
      load([dirname,list.mat{idx}],'numFrameErr','numFramesSimulated','EbNodB_vec')
      numFrameErr_cum         = numFrameErr_cum + numFrameErr;
      numFramesSimulated_cum  = numFramesSimulated_cum + numFramesSimulated;
    end
  end
  if (fileCount ~= fileCountExpected)
    keyboard;
  end
  if any(numFrameErr_cum > 0)
    %disp(numFrameErr_cum)
    semilogy(EbNodB_vec,numFrameErr_cum./numFramesSimulated_cum,plotstyle{bsn_idx},'DisplayName',basefilename{bsn_idx})
    hold on, grid on
  end
end

xlabel('Eb/No [dB]')
ylabel('FER')
