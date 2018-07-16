function compare_outputs(alg)
base_dir = 'F:/Research/SinusProject/Seth_code';

% system(sprintf(...
%     '%s/cisstICP/cissticp/cisstICP_App/build/Release/ICP_App --alg DIMLP --out testingforrelease', ...
%         base_dir));

parent_dir = sprintf('%s/cisstICP/cissticp/test_data/LastRun_%s', base_dir, alg);
original = sprintf('%s/testing/', parent_dir);
modified = sprintf('%s/testingforrelease', parent_dir);
o_res = read_outputs(sprintf('%s/SaveIterations.txt', original));
m_res = read_outputs(sprintf('%s/SaveIterations.txt', modified));
if norm(o_res - m_res) < 0.1
    disp('pass!')
else
    disp(norm(o_res - m_res));
    disp('fail!')
end

fclose('all');
clear all

function [res] = read_outputs(filename)

    modes = 3;
    fid=fopen(filename);
    N=textscan(fid, '%s', 'delimiter', '\n');
    sN = size(N{1,1},1);
    str = N{1,1}{sN-3};
    C=strsplit(str);
    res(1) = str2double(C(4));          % rotation
    res(2) = str2double(C(6));          % translation
    if size(C,2) > 6
        for s = 1:modes
            res(2+s) = str2double(C(7+s));  % shape parameters 
        end
    end
    str = N{1,1}{sN-1};
    C=strsplit(str,{' ','-',':',')'});
    res(6) = str2double(C(5));          % mean residual error
    res(7) = str2double(C(7));          % standard deviation
    
end
end