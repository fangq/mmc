
% full field illumination pattern
pat0 = ones(40);

% pattern with left half bright
pat1 = zeros(40);
pat1(1:20,:) = 1;

% pattern with top half bright
pat2 = zeros(40);
pat2(:,1:20) = 1;

% pattern with bright square in the middle
pat3 = zeros(40);
pat3(11:30,11:30) = 1;

% pattern with a bright cross
pat4 = zeros(40);
pat4(16:25,:) = 1;
pat4(:,16:25) = 1;

fid=fopen('patterns.pat','wb','ieee-le');
fwrite(fid,pat0','float');
fwrite(fid,pat1','float');
fwrite(fid,pat2','float');
fwrite(fid,pat3','float');
fwrite(fid,pat4','float');
fclose(fid);