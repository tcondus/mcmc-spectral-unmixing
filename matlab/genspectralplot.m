% Edits by (tc): Adjustments to the plots based on personal preference.
% For combined S and L data (i.e., ~0.4 to 2.5 microns), the plots are
% divided for better clarity. However, the bands are currently hard coded.
function [R2, rmse, lt, wt, wmi, wifi, f] = genspectralplot(ab, d, ind, all_ab, all_d, all_ind, R_vs_SSA, D)

set(0, 'defaultAxesFontSize', 18)
set(0, 'defaultTextFontSize', 18)
set(0, 'defaultfigurecolor', [1 1 1])

ABm = ab;
Dm = d;
INDm = ind;

lt = D.lam_SPEC;
wt = D.R_SPEC;

[Mm, f] = gen_synthetic(ABm, Dm, INDm, R_vs_SSA, D);
lm = Mm(:, 1);
wm = Mm(:, 2);
wif = Mm(:, 3:end);

wmi = interp1(lm, wm, lt);
wifi = interp1(lm, wif, lt);
yb = 1/length(wmi)*sum(wt);
inddd = find(~isnan(wmi));
wmi = wmi(inddd);
wifi = wifi(inddd, 1:end);
wt = wt(inddd);
lt = lt(inddd);
SStot = sum((wt-yb).^2);
SSres = sum((wt-wmi).^2);
R2 = 1 - SSres/SStot;
[T, ~] = size(inddd);
rmse = sqrt(SSres / T);

figure
if round(max(lt)*1e6) == 1 || round(min(lt)*1e6) == 1 % If S data or L data only
    subplot(2, 2, 1, 'align')
else
    subplot(2, 3, 1, 'align')
end
hold on

% Plot all model spectra at mixed spectrum resolution.
[~, all_ab_cols] = size(all_ab);
for i = 1:all_ab_cols
    all_ABm = all_ab(:, i);
    all_Dm = all_d(:, i);
    all_INDm = all_ind(:, i);
    
    % Shorten the vectors when num_endmembers < max_endmembers.
    all_ABm = all_ABm(all_ABm ~= 0);
    all_Dm = all_Dm(all_Dm ~= 0);
    all_INDm = all_INDm(all_INDm ~= 0);
    
    all_Mm = gen_synthetic(all_ABm, all_Dm, all_INDm, R_vs_SSA, D);
    all_lm = all_Mm(:, 1);
    all_wm = all_Mm(:, 2);
    all_wmi = interp1(all_lm, all_wm, lt);
    if ~(round(max(lt)*1e6) == 1 || round(min(lt)*1e6) == 1) % If whole spectrum (from 0.4 to 2.5 microns)
        plot(lt(1:79).*1e6, all_wmi(1:79), 'Color', [0.75 0.75 0.75], 'LineWidth', 1)
    else
        plot(lt.*1e6, all_wmi, 'Color', [0.75 0.75 0.75], 'LineWidth', 1)
    end
    %plot(lt.*1e6, all_wmi, 'k', 'LineWidth', 1)
end
if ~(round(max(lt)*1e6) == 1 || round(min(lt)*1e6) == 1) % If whole spectrum (from 0.4 to 2.5 microns)
    p1 = plot(lt(1:79).*1e6, wt(1:79), 'r', 'LineWidth', 3);
    p2 = plot(lt(1:79).*1e6, wmi(1:79), 'k', 'LineWidth', 2); % Plot model spectrum at mixed spectrum resolution.
else
    p1 = plot(lt.*1e6, wt, 'r', 'LineWidth', 3);
    p2 = plot(lt.*1e6, wmi, 'k', 'LineWidth', 2); % Plot model spectrum at mixed spectrum resolution.
end
% wc = zeros(length(wmi), 1);
% for i = 1:length(ABm)
%     wc = wc + (ABm(i) * wifi(:, i));
% end
% p3 = plot(lt.*1e6, wc, 'b', 'LineWidth', 2); % Plot checkerboard spectrum at mixed spectrum resolution.
%plot(lm.*1e6, wm, 'k', 'LineWidth', 2) % Plot model spectrum at n and k resolution.
legend([p1 p2], sprintf('data R^2 = %.4f', R2), sprintf('model RMSE = %.4f', rmse), 'Location', 'NorthOutside')
xlabel('Wavelength, \lambda (\mum)')
if R_vs_SSA
    ylabel('Reflectance')
else
    ylabel('Single Scattering Albedo')
end
if ~(round(max(lt)*1e6) == 1 || round(min(lt)*1e6) == 1) % If whole spectrum (from 0.4 to 2.5 microns)
    xlim([min(lt(1:79).*1e6) max(lt(1:79).*1e6)])
    ylim([min(wt(1:79))-0.01 max(wt(1:79))+0.01])
else
    xlim([min(lt.*1e6) max(lt.*1e6)])
    ylim([min(wt)-0.01 max(wt)+0.01])
end
%xlim([min([min(lm.*1e6) min(lt.*1e6)]) max([max(lm.*1e6) max(lt.*1e6)])])
%ylim([-inf inf])

set(gcf, 'color', 'w');
if round(max(lt)*1e6) == 1 || round(min(lt)*1e6) == 1 % If S data or L data only
    subplot(2, 2, 3, 'align')
else
    subplot(2, 3, 4, 'align')
end
if ~(round(max(lt)*1e6) == 1 || round(min(lt)*1e6) == 1) % If whole spectrum (from 0.4 to 2.5 microns)
    plot(lt(1:79).*1e6,(wt(1:79)./wmi(1:79)-1).*100, 'r', 'LineWidth', 2)
    hold on
    plot(lt(1:79).*1e6, zeros(size(lt(1:79))), 'k', 'LineWidth', 1)
else
    plot(lt.*1e6,(wt./wmi-1).*100, 'r', 'LineWidth', 2)
    hold on
    plot(lt.*1e6, zeros(size(lt)), 'k', 'LineWidth', 1)
end
xlabel('Wavelength, \lambda (\mum)')
ylabel('Residual (%)')
if ~(round(max(lt)*1e6) == 1 || round(min(lt)*1e6) == 1) % If whole spectrum (from 0.4 to 2.5 microns)
    xlim([min(lt(1:79).*1e6) max(lt(1:79).*1e6)])
    if abs(100.*(wt(1:79)./wmi(1:79)-1)) < 5
        ylim([-5 5])
    end
else
    xlim([min(lt.*1e6) max(lt.*1e6)])
    if abs(100.*(wt./wmi-1)) < 5
        ylim([-5 5])
    end
end
%xlim([min([min(lm.*1e6) min(lt.*1e6)]) max([max(lm.*1e6) max(lt.*1e6)])])

if ~(round(max(lt)*1e6) == 1 || round(min(lt)*1e6) == 1) % If whole spectrum (from 0.4 to 2.5 microns)
    subplot(2, 3, 2, 'align')
    hold on
    % Plot all model spectra at mixed spectrum resolution.
    [~, all_ab_cols] = size(all_ab);
    for i = 1:all_ab_cols
        all_ABm = all_ab(:, i);
        all_Dm = all_d(:, i);
        all_INDm = all_ind(:, i);

        % Shorten the vectors when num_endmembers < max_endmembers.
        all_ABm = all_ABm(all_ABm ~= 0);
        all_Dm = all_Dm(all_Dm ~= 0);
        all_INDm = all_INDm(all_INDm ~= 0);

        all_Mm = gen_synthetic(all_ABm, all_Dm, all_INDm, R_vs_SSA, D);
        all_lm = all_Mm(:, 1);
        all_wm = all_Mm(:, 2);
        all_wmi = interp1(all_lm, all_wm, lt);
        plot(lt(79:end).*1e6, all_wmi(79:end), 'Color', [0.75 0.75 0.75], 'LineWidth', 1)
        %plot(lt.*1e6, all_wmi, 'k', 'LineWidth', 1)
    end
    p1 = plot(lt(79:end).*1e6, wt(79:end), 'r', 'LineWidth', 3);
    p2 = plot(lt(79:end).*1e6, wmi(79:end), 'k', 'LineWidth', 2); % Plot model spectrum at mixed spectrum resolution.

    % wc = zeros(length(wmi), 1);
    % for i = 1:length(ABm)
    %     wc = wc + (ABm(i) * wifi(:, i));
    % end
    % p3 = plot(lt.*1e6, wc, 'b', 'LineWidth', 2); % Plot checkerboard spectrum at mixed spectrum resolution.
    %plot(lm.*1e6, wm, 'k', 'LineWidth', 2) % Plot model spectrum at n and k resolution.
    legend([p1 p2], sprintf('data R^2 = %.4f', R2), sprintf('model RMSE = %.4f', rmse), 'Location', 'NorthOutside')
    xlabel('Wavelength, \lambda (\mum)')
    if R_vs_SSA
        ylabel('Reflectance')
    else
        ylabel('Single Scattering Albedo')
    end
    xlim([min(lt(79:end).*1e6) max(lt(79:end).*1e6)])
    ylim([min(wt(79:end))-0.01 max(wt(79:end))+0.01])
    %xlim([min([min(lm.*1e6) min(lt.*1e6)]) max([max(lm.*1e6) max(lt.*1e6)])])
    %ylim([-inf inf])
    
    set(gcf, 'color', 'w');
    subplot(2, 3, 5, 'align')
    plot(lt(79:end).*1e6,(wt(79:end)./wmi(79:end)-1).*100, 'r', 'LineWidth', 2)
    hold on
    plot(lt(79:end).*1e6, zeros(size(lt(79:end))), 'k', 'LineWidth', 1)
    xlabel('Wavelength, \lambda (\mum)')
    ylabel('Residual (%)')
    xlim([min(lt(79:end).*1e6) max(lt(79:end).*1e6)])
    %xlim([min([min(lm.*1e6) min(lt.*1e6)]) max([max(lm.*1e6) max(lt.*1e6)])])
    if abs(100.*(wt(79:end)./wmi(79:end)-1)) < 5
        ylim([-5 5])
    end
else

if round(max(lt)*1e6) == 1 || round(min(lt)*1e6) == 1 % If S data or L data only
    subplot(2, 2, 2, 'align')
else
    subplot(2, 3, 3, 'align')
end
hold on
[ABm_rows, ~] = size(ABm);
for i = 1:ABm_rows
    b = barh(ABm_rows-(i-1), ABm(i));
    xtips1 = b(1).YData + 0.01;
    ytips1 = b(1).XData;
    labels1 = strcat(string(round(b(1).YData*100, 2)), '%');
    text(xtips1, ytips1, labels1, 'VerticalAlignment', 'middle')
end
yticklabels(repmat({' '}, 1, ABm_rows));

N = length(ind);
labels = strings(1, 5);
for i=1:N
    index = ind(i);
    labels(1, i) = D.labels(index);
end

hFig = gcf;
if round(max(lt)*1e6) == 1 || round(min(lt)*1e6) == 1 % If S data or L data only
    %h = legend(D.labels);
    h = legend(labels);
    rect = [0.60, 0.15, 0.25, 0.25]; % First number: position from left. Second number: position from bottom.
    set(h, 'Position', rect)
    
    hFig.Units = 'inches';
    set(hFig, 'Position', [1 1 12 8])
else
    set(hFig, 'Position', [100 100 500 1000])
end
set(hFig, 'Visible', 'off'); % Prevent MATLAB from opening the figure in a new window.

end
