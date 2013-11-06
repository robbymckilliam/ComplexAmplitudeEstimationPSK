	if debugplot
    subplot(3,1,1);
    cla;
    for i=1:length(SoftDec)
      if SoftDec(i).Acquired
        if SoftDec(i).finish
          plot(SoftDec(i).Fo*ModemParams.Fs/1e3, SoftDec(i).Td/ModemParams.Fs, 'go', 'MarkerSize', 10, 'LineWidth', 2);
        elseif SoftDec(i).success
          plot(SoftDec(i).Fo*ModemParams.Fs/1e3, SoftDec(i).Td/ModemParams.Fs, 'gx', 'MarkerSize', 10, 'LineWidth', 2);
        else
          plot(SoftDec(i).Fo*ModemParams.Fs/1e3, SoftDec(i).Td/ModemParams.Fs, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
        end
        text(SoftDec(i).Fo*ModemParams.Fs/1e3, SoftDec(i).Td/ModemParams.Fs + 0.5/ModemParams.SymbolRate,...
             sprintf('   %.2f dB', 10*log10(SoftDec(i).SigPowEst/SoftDec(i).NoiseVar)), 'Rotation', 90, 'FontSize', 10);
        text(SoftDec(i).Fo*ModemParams.Fs/1e3, SoftDec(i).Td/ModemParams.Fs - 0.5/ModemParams.SymbolRate,...
             sprintf('%.i', i), 'FontSize', 10, 'HorizontalAlignment', 'center');
      end
      hold on;
    end
    axis([-ModemParams.Fs/2/1e3 ModemParams.Fs/2/1e3 0 max([1.2*ModemParams.MaxTd 1e-3]) ]);
    title(sprintf('Iteration %i\n%i acquired, %i decoded, %i finished', m, sum([SoftDec.Acquired]), sum([SoftDec.success]), sum([SoftDec.finish])));
    xlabel('Frequency (kHz)');
    ylabel('Time (s)');
    
    yrange = [-80 -20];
    subplot(3,1,2);
    cla;
    hold off;
    psd(spectrum.periodogram, NoiseHypothesis, 'Fs', ModemParams.Fs, 'CenterDC', true);
    ylim(yrange);
%     title(sprintf('%.2f dB', 10*log10(mean(abs(NoiseHypothesis.^2))/ModemParams.Fs)));
    title('cancelled signal');
    hold on;
%     mask  = fftshift(STAcqParams.FFTMask);
    if isfield(STAcqParams, 'FFTMask')
      mask  = STAcqParams.FFTMask;
      mi    = 1;
      while(mi<length(mask))
        % find exclusion start
        idx1 = find(1-mask(mi:end), 1)+mi-1;
        if ~isempty(idx1)
          idx2 = find(mask(idx1:end), 1)+idx1-1;
          if isempty(idx2)
            idx2 = length(mask);
          end
          patch(([idx1 idx1 idx2 idx2]-length(mask)/2)*ModemParams.Fs/length(mask)/1e3, [yrange fliplr(yrange)], [0 0 0 0], 'facecolor', 'r', 'linestyle', 'none', 'facealpha', 0.5);
        else
          break;
        end
        mi = idx2 + 1;
      end
    end
    
    if m==1
      subplot(3,1,3);
      cla;
      psd(spectrum.periodogram, RxSamples, 'Fs', ModemParams.Fs, 'CenterDC', true);
      ylim(yrange);
%       title(sprintf('%.2f dB', 10*log10(mean(abs(RxSamples.^2))/ModemParams.Fs)));
      title('received signal');
    end
    
    drawnow;
    
%     set(gcf,'PaperUnits','centimeters','PaperSize',[29,20],'PaperPosition',[0 0 29 20]);
%     print('-dpng', '-r300', sprintf('st1receiver%03i.png', m));
  end
  