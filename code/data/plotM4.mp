prologues := 2;
filenametemplate "%j-%c.mps";

input TEX; 
TEXPRE("%&latex" & char(10) & 
"\documentclass{amsart}\begin{document}"); 
TEXPOST("\end{document}"); 

input graph;
input sarith;
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

ls:=2;
label.rt(btex $|P|=L/32$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw circle shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=3;
label.rt(btex $|P|=0$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw cross shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=4;
label.rt(btex $|P|=L$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw times shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=5;
label.rt(btex Viterbi Viterbi etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups) dashed evenly;

ls:=6;
label.rt(btex Theory etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);

ls:=7;
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
label.rt(btex $|P|=L/8$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw dot shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

%ls:=2;
%label.rt(btex $|P|=L/32$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
%draw circle shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

%ls:=3;
%label.rt(btex $|P|=0$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
%draw cross shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=1;
label.rt(btex $|P|=L$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw times shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=2;
label.rt(btex Theory etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);

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
  glabel.lft(btex variance $\hat{\rho}$ etex, OUT) rotated 90 shifted (1.5cm,1cm);
  
  %setrange(whatever,whatever, whatever,5);
  %setrange(4.9,"1e-13", 22,2);

%  gdraw (0, "1e-13")--(0, 2) dashed evenly; %glabel.bot(btex AWGN channel capacity etex, 0) rotated 90;

  gdraw "cltM4L2048absP2048a"; gdraw "MackenthunCoherentM4L2048absP2048aunb" plot(times); glabel.llft(btex $L = 2048$ etex, 1000);
  gdraw "cltM4L2048absP1024a"; gdraw "MackenthunCoherentM4L2048absP1024aunb" plot(square);
%  gdraw "cltM4L2048absP512a"; gdraw "MackenthunCoherentM4L2048absP512aunb" plot(dot);
  %gdraw "cltM4L2048absP128a"; gdraw "MackenthunCoherentM4L2048absP128aunb" plot(circle);
  gdraw "cltM4L2048absP0a"; gdraw "MackenthunNonCoherentM4L2048aunb" plot(cross);

  gdraw "cltM4L256absP256a"; gdraw "MackenthunCoherentM4L256absP256aunb" plot(times); glabel.llft(btex $L = 256$ etex, 1000);
  gdraw "cltM4L256absP128a"; gdraw "MackenthunCoherentM4L256absP128aunb" plot(square);
  %gdraw "cltM4L256absP512a"; gdraw "MackenthunCoherentM4L256absP512aunb" plot(dot);
  %gdraw "cltM4L256absP128a"; gdraw "MackenthunCoherentM4L256absP128aunb" plot(circle);
  gdraw "cltM4L256absP0a"; gdraw "MackenthunNonCoherentM4L256aunb" plot(cross);

  gdraw "cltM4L32absP32a"; gdraw "MackenthunCoherentM4L32absP32aunb" plot(times); glabel.llft(btex $L = 32$ etex, 1000);
  gdraw "cltM4L32absP16a"; gdraw "MackenthunCoherentM4L32absP16aunb" plot(square);
 % gdraw "cltM4L32absP512a"; gdraw "MackenthunCoherentM4L32absP512aunb" plot(dot);
  %gdraw "cltM4L32absP128a"; gdraw "MackenthunCoherentM4L32absP128aunb" plot(circle);
  gdraw "cltM4L32absP0a"; gdraw "MackenthunNonCoherentM4L32aunb" plot(cross);
  
draw legendamp shifted(legendx,legendy);

endgraph;

%clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight)--(plotwidth+5cm,plotheight)--(plotwidth+5cm,-10cm)--cycle

endfig;


%%%%%%%%%%%%
%%%% Begin plot of phase error

beginfig(2);
draw begingraph(plotwidth, plotheight);
 setcoords(linear,log);
 glabel.bot(btex SNR (dB) etex, OUT) shifted (0cm,1.3cm);
  glabel.lft(btex MSE $\hat{\theta}$ etex, OUT) rotated 90 shifted (1.5cm,1cm);
  
  setrange(whatever,"9e-7", whatever,0.99);
  %setrange(4.9,"1e-13", 22,2);

  gdraw "cltM4L4096absP4096p"; gdraw "MackenthunCoherentM4L4096absP4096p" plot(times); %glabel.llft(btex $L = 4096$ etex, 1000);
  gdraw "cltM4L4096absP2048p"; gdraw "MackenthunCoherentM4L4096absP2048p" plot(square);
  gdraw "cltM4L4096absP512p"; gdraw "MackenthunCoherentM4L4096absP512p" plot(dot);
  gdraw "cltM4L4096absP128p"; gdraw "MackenthunCoherentM4L4096absP128p" plot(circle);
  gdraw "cltM4L4096absP0p"; gdraw "MackenthunNonCoherentM4L4096p" plot(cross);
  gdraw "ViterbiViterbiM4L4096p" dashed evenly;
  gdraw "QPSKCRB4096" dashed dd;

 % gdraw "cltM4L2048absP2048p"; gdraw "MackenthunCoherentM4L2048absP2048p" plot(times); %glabel.llft(btex $L = 2048$ etex, 1000);
 % gdraw "cltM4L2048absP1024p"; gdraw "MackenthunCoherentM4L2048absP1024p" plot(dot);
 % gdraw "cltM4L2048absP256p"; gdraw "MackenthunCoherentM4L2048absP256p" plot(square);
%  gdraw "cltM4L2048absP64p"; gdraw "MackenthunCoherentM4L2048absP64p" plot(cross);
%  gdraw "cltM4L2048absP0p"; gdraw "MackenthunNonCoherentM4L2048p" plot(cross);

  %gdraw "cltM4L128absP128p"; gdraw "MackenthunCoherentM4L128absP128p" plot(times); glabel.llft(btex $L = 128$ etex, 1000);
  %gdraw "cltM4L128absP64p"; gdraw "MackenthunCoherentM4L128absP64p" plot(dot);
  %gdraw "cltM4L128absP16p"; gdraw "MackenthunCoherentM4L128absP16p" plot(square);
  %gdraw "cltM4L128absP4p"; gdraw "MackenthunCoherentM4L128absP4p" plot(cross);
  %gdraw "cltM4L128absP0p"; gdraw "MackenthunNonCoherentM4L128p" plot(cross);

  
draw legend shifted(legendx,legendy);

endgraph;

%clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight)--(plotwidth+5cm,plotheight)--(plotwidth+5cm,-10cm)--cycle

endfig;



%%%%%%%%%%%%
%%%% Begin plot of phase error multiple L

beginfig(3);
draw begingraph(plotwidth, plotheight);
 setcoords(linear,log);
 glabel.bot(btex SNR (dB) etex, OUT) shifted (0cm,1.3cm);
  glabel.lft(btex MSE $\hat{\theta}$ etex, OUT) rotated 90 shifted (1.5cm,1cm);
  
  setrange(whatever,"9e-7", whatever,5);
  %setrange(4.9,"1e-13", 22,2);

%  gdraw "cltM4L4096absP4096p"; gdraw "MackenthunCoherentM4L4096absP4096p" plot(times); %glabel.llft(btex $L = 4096$ etex, 1000);
%  gdraw "cltM4L4096absP2048p"; gdraw "MackenthunCoherentM4L4096absP2048p" plot(dot);
%  gdraw "cltM4L4096absP512p"; gdraw "MackenthunCoherentM4L4096absP512p" plot(square);
%  gdraw "cltM4L4096absP128p"; gdraw "MackenthunCoherentM4L4096absP128p" plot(circle);
%  gdraw "cltM4L4096absP0p"; gdraw "MackenthunNonCoherentM4L4096p" plot(cross);

  gdraw "cltM4L2048absP2048p"; gdraw "MackenthunCoherentM4L2048absP2048p" plot(times); label.llft(btex $L = 2048$ etex,(0,0)) rotated -30 shifted (9cm,0.6cm);
%  gdraw "cltM4L2048absP1024p"; gdraw "MackenthunCoherentM4L2048absP1024p" plot(square);
  gdraw "cltM4L2048absP256p"; gdraw "MackenthunCoherentM4L2048absP256p" plot(dot);
%  gdraw "cltM4L2048absP64p"; gdraw "MackenthunCoherentM4L2048absP64p" plot(circle);
%  gdraw "cltM4L2048absP0p"; gdraw "MackenthunNonCoherentM4L2048p" plot(cross);

  gdraw "cltM4L256absP256p"; gdraw "MackenthunCoherentM4L256absP256p" plot(times); label.llft(btex $L = 256$ etex,(0,0)) rotated -30 shifted (9cm,1.7cm);
%  gdraw "cltM4L256absP128p"; gdraw "MackenthunCoherentM4L256absP128p" plot(square);
  gdraw "cltM4L256absP32p"; gdraw "MackenthunCoherentM4L256absP32p" plot(dot);
%  gdraw "cltM4L256absP64p"; gdraw "MackenthunCoherentM4L256absP64p" plot(circle);
%  gdraw "cltM4L256absP0p"; gdraw "MackenthunNonCoherentM4L256p" plot(cross);

  gdraw "cltM4L32absP32p"; gdraw "MackenthunCoherentM4L32absP32p" plot(times); label.llft(btex $L = 32$ etex,(0,0)) rotated -30 shifted (9cm,2.8cm);
%  gdraw "cltM4L32absP128p"; gdraw "MackenthunCoherentM4L32absP128p" plot(square);
  gdraw "cltM4L32absP4p"; gdraw "MackenthunCoherentM4L32absP4p" plot(dot);
%  gdraw "cltM4L32absP64p"; gdraw "MackenthunCoherentM4L32absP64p" plot(circle);
%  gdraw "cltM4L32absP0p"; gdraw "MackenthunNonCoherentM4L32p" plot(cross);

  
draw legendphcomp shifted(legendx,legendy);

endgraph;

%clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight)--(plotwidth+5cm,plotheight)--(plotwidth+5cm,-10cm)--cycle

endfig;


end;