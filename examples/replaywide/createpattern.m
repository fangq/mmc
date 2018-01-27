
% full field illumination pattern
pat0 = ones(40);

% pattern with left half dark
pat1 = pat0;
pat1(1:20,:) = 0;

% pattern with top half dark
pat2 = pat0;
pat2(:,1:20) = 0;

fid=fopen('pattern0.pat','wb','ieee-le');
fwrite(fid,pat0','float');
fclose(fid);

fid=fopen('pattern1.pat','wb','ieee-le');
fwrite(fid,pat1','float');
fclose(fid);

fid=fopen('pattern2.pat','wb','ieee-le');
fwrite(fid,pat2','float');
fclose(fid);