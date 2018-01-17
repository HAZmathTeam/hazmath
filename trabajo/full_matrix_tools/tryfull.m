 m=4;n=11;p=7; a=randn(m,n);b=randn(n,p);
xn=randn(n,1);xm=randn(m,1); ym=zeros(m,1);yn=zeros(n,1); c=zeros(m,p);
fid=fopen('abc.input','w');
fprintf(fid,'%23.16e ',a');
 fprintf(fid,'\n');
 fprintf(fid,'%23.16e ',b');
fprintf(fid,'\n');
fprintf(fid,'%23.16e ',c');
 fprintf(fid,'\n');
 fprintf(fid,'%23.16e ',xm);
 fprintf(fid,'\n');
 fprintf(fid,'%23.16e ',xn);
 fprintf(fid,'\n');
 fprintf(fid,'%23.16e ',ym);
 fprintf(fid,'\n');
 fprintf(fid,'%23.16e ',yn);
 fclose(fid)
