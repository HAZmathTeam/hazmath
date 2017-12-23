function a = plotgraf2(x,y,nfile,nline)
xy = [x,y];


             fortk = ['fort.' int2str(nfile)];
 	     filename = [fortk]
             if ~exist(filename), clear, close all, break, end
             eval(['load ' filename])
a = sparse(fort(1,:)',fort(2,:)',fort(3,:)');
clear fort
%a = a + a';
if(nline==1)
gplot(a,xy,'k');
else if(nline==0)
gplot(a,xy,'k:');
else
gplot(a,xy,'c:');
end
end
return