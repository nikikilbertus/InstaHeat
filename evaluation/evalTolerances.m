%% evaluate accuracy and precision for runs with different tolerances

%% setup (user input)
% setup different tolerances
apot = [4 5 6 7 8 9 10 11 12 13 14];
rpot = [4 5 6 7 8 9 10 11 12 13 14];
% construct the file names: prefix, suffix, indexset
pre = '';
suf = '';
%figure offset
os = 0;

%% run
loadDsets;

reltols = 10.^(-rpot);
abstols = 10.^(-apot);

getname = @(x,y) [pre num2str(x) '_' num2str(y) suf '.h5'];
nr = length(rpot); na = length(apot); nn = nr * na;

runtime = zeros([nr, na]);
steps = zeros([nr, na]);
phimeanerrl2 = zeros([nr, na]);
phistderrl2 = zeros([nr, na]);
rhormserrl2 = zeros([nr, na]);
as = zeros([nr, na]);
cstrl2 = zeros([nr, na]);

clear dsets
require('a');
name = getname(min(rtol),min(atol));
readDsets;
ashare = a;

clear dsets
require('phiS', 'a', 'constraints', 'rhoS')
name = getname(max(rpot),max(apot));
readDsets;
aref = a;
phimeanref = spline(a,phimean,ashare);
phistdref = spline(a,phistd,ashare);
rhormsref = spline(a,rhorms,ashare);
cstrl2ref = spline(a,constraints(:,1),ashare);

clear dsets
require('phiS','rhoS','constraints','steps_total','tolerances','runtime_total');
for i = 1:nr
    for j = 1:na
        name = getname(rpot(i),apot(j));
        readDsets;
        if reltols(i) ~= tols(1) || abstols(j) ~= tols(2)
            error('didnt load the correct file')
        end
        runtime(i,j) = runtime_total;
        steps(i,j) = steps_total;
        as(i,j) = -log10(abs((a(end) - aref(end))/aref(end)));
        
        phimean = spline(a,phimean,ashare);
        phistd = spline(a,phistd,ashare);
        rhorms = spline(a,rhorms,ashare);
        cstr = spline(a,constraints(:,1),ashare);
        
        phimeanerrl2(i,j) = norm(phimeanref - phimean);
        phistderrl2(i,j) = norm(phistdref - phistd);
        rhormserrl2(i,j) = norm(rhormsref - rhorms);
        cstrl2(i,j) = -log10(abs((cstr - cstrl2ref) / cstrl2ref));
 
%         I = (a>0.9*aref(end));
%         Iref = (aref>0.9*aref(end));
%         plot(a(I),phi(1,I).*a(I)'.^(3/2),aref(Iref),phiref(Iref).*aref(Iref)'.^(3/2)); shg; pause;
    end
end
% semilogx(reltols, time, 'linewidth',2); xlabel('rel tol'); ylabel('time [s]');
% legend('abs: 1e-6', 'abs: 1e-8', 'abs: 1e-10', 'abs: 1e-12'); shg;
% figure
% semilogx(reltols, steps, 'linewidth',2); xlabel('rel tol'); ylabel('#steps');
% legend('abs: 1e-6', 'abs: 1e-8', 'abs: 1e-10', 'abs: 1e-12'); shg;
% figure
% semilogx(abstols, time', 'linewidth',2); xlabel('abs tol'); ylabel('time [s]');
% legend('rel: 1e-6', 'rel: 1e-8', 'rel: 1e-10', 'rel: 1e-12'); shg;
% figure
% semilogx(abstols, steps', 'linewidth',2); xlabel('abs tol'); ylabel('#steps');
% legend('rel: 1e-6', 'rel: 1e-8', 'rel: 1e-10', 'rel: 1e-12'); shg;

figure
bar3(as); set(gca,'XTickLabel',abstols); set(gca,'YTickLabel',reltols);
xlabel('atol'); ylabel('rtol'); zlabel('-log10 (a_{f} - a_{f}^{ref}) / a_{f}^{ref}');
figure
bar3(steps); set(gca,'XTickLabel',abstols); set(gca,'YTickLabel',reltols);
xlabel('atol'); ylabel('rtol'); zlabel('#steps');
figure
bar3(cstrl2); set(gca,'XTickLabel',abstols); set(gca,'YTickLabel',reltols);
xlabel('atol'); ylabel('rtol'); zlabel('hamiltonian constraint norm');
figure
bar3(phimeanerrl2); set(gca,'XTickLabel',abstols); set(gca,'YTickLabel',reltols);
xlabel('atol'); ylabel('rtol'); zlabel('-log10 <\phi> error l_{2}');
figure
bar3(phistderrl2); set(gca,'XTickLabel',abstols); set(gca,'YTickLabel',reltols);
xlabel('atol'); ylabel('rtol'); zlabel('-log10 std \phi error l_{2}');
figure
bar3(rhormserrl2); set(gca,'XTickLabel',abstols); set(gca,'YTickLabel',reltols);
xlabel('atol'); ylabel('rtol'); zlabel('-log10 std \phi error l_{2}');