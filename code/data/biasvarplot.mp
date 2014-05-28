prologues := 2;
filenametemplate "%j-%c.mps";

input graph;
input sarith;

input TEX; 
TEXPRE("%&latex" & char(10) & 
"\documentclass{minimal}\begin{document}"); 
TEXPOST("\end{document}"); 

%verbatimtex
%&latex
%\documentclass[8pt]{extarticle}
%  \begin{document}
%  etex

%init_numbers(btex $-$ etex,btex $1$ etex, btex ${\times}10$ etex, btex ${}^-$ etex, btex ${}^2$ etex);
  
picture ORIGINAL_PICTURE;
ORIGINAL_PICTURE := currentpicture;

picture dd; dd :=
dashpattern(on 4bp off 2bp on 1bp off 2bp);

picture diamond;
draw unitsquare shifted(-0.5,-0.5) scaled 0.17cm rotated 45;
diamond := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture square;
draw unitsquare shifted(-0.5,-0.5) scaled 0.15cm;
square := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture filledsquare;
fill unitsquare shifted(-0.5,-0.5) scaled 0.12cm;
filledsquare := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture circle; 
draw fullcircle scaled 0.19cm;
circle := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture cross; 
draw (-0.09cm,0.0)--(0.09cm,0.0);
draw (0,-0.09cm)--(0,0.09cm);
cross := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture dot; 
fill fullcircle scaled 0.15cm;
dot := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture times; times := cross rotated 45;


%%%%%%%%%%%%%
%%%%%Legend picture

picture legend;
sx = 0.2cm;
sh = 0.2cm;
ups = 0.4cm;

ls:=0;
label.rt(btex $|P|=L/2$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw square shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=1;
label.rt(btex $|P|=L/8$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw dot shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

%ls:=2;
%label.rt(btex $|P|=L/32$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
%draw circle shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=2;
label.rt(btex $|P|=0$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw cross shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=3;
label.rt(btex $|P|=L$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw times shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=4;
label.rt(btex Viterbi Viterbi etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups) dashed evenly;

ls:=5;
label.rt(btex Theory etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);

ls:=6;
label.rt(btex CRB etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups) dashed dd;

legend := currentpicture; currentpicture := ORIGINAL_PICTURE;


picture legendamp;

ls:=0;
label.rt(btex $|P|=L/2$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw square shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

%ls:=1;
%label.rt(btex $|P|=L/8$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
%draw dot shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

%ls:=2;
%label.rt(btex $|P|=L/32$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
%draw circle shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=1;
label.rt(btex $|P|=0$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw cross shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=2;
label.rt(btex $|P|=L$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw times shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=3;
label.rt(btex Theory etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);


legendamp := currentpicture; currentpicture := ORIGINAL_PICTURE;


picture legendphcomp;

%ls:=0;
%label.rt(btex $|P|=L/2$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
%draw square shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=0;
label.rt(btex bias squared etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
%draw dot shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

%ls:=2;
%label.rt(btex $|P|=L/32$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
%draw circle shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

%ls:=3;
%label.rt(btex $|P|=0$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
%draw cross shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=1;
label.rt(btex variance etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups) dashed evenly;
%draw times shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

%ls:=2;
%label.rt(btex Theory etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);

legendphcomp := currentpicture; currentpicture := ORIGINAL_PICTURE;

plotheight = 8cm;
plotwidth = 10cm;

legendx = 0.1cm;
legendy = 0.1cm;


%%%%%%%%%%%%
%%%% Begin plot of amplitude error

beginfig(1);
draw begingraph(plotwidth, plotheight);
 setcoords(linear,log);
 glabel.bot(btex SNR (dB) etex, OUT) shifted (0cm,1.3cm);
  %glabel.lft(btex MSE $\hat{\rho}$ etex, OUT) rotated 90 shifted (1.5cm,0cm);

  
  setrange(-16,"1e-9", 12, 50); 
  %setrange(4.9,"1e-13", 22,2);

%  gdraw (0, "1e-13")--(0, 2) dashed evenly; %glabel.bot(btex AWGN channel capacity etex, 0) rotated 90;

  %gdraw "cltM4L4096absP4096b"; gdraw "cltM4L4096absP4096a"; 
  gdraw "cltM4L4096absP2048b"; glabel.bot(btex $|P|=L/2$ etex, 20);
  gdraw "cltM4L4096absP2048a" dashed evenly; glabel.top(btex $|P|=L/2$ etex, 20);
  
  gdraw "cltM4L4096absP0b"; glabel.top(btex $|P|=0$ etex, 20);
  gdraw "cltM4L4096absP0a" dashed evenly; glabel.bot(btex $|P|=0$ etex, 20);
  
draw legendphcomp shifted(legendx,legendy);

endgraph;

%clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight)--(plotwidth+5cm,plotheight)--(plotwidth+5cm,-10cm)--cycle

endfig;


end;
