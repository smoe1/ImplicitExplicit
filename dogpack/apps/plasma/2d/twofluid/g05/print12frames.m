function print12frames(n)

sn = num2str(n);

% File names
fname01 = 'rhoi_xx_color.jpg';
if n<10
 fname01(6)='0';
 fname01(7)=sn(1);
else
 fname01(6)=sn(1);
 fname01(7)=sn(2);
end

fname02 = 'rhoi_xx_contour.jpg';
if n<10
 fname02(6)='0';
 fname02(7)=sn(1);
else
 fname02(6)=sn(1);
 fname02(7)=sn(2);
end

fname03 = 'rhoe_xx_color.jpg';
if n<10
 fname03(6)='0';
 fname03(7)=sn(1);
else
 fname03(6)=sn(1);
 fname03(7)=sn(2);
end

fname04 = 'rhoe_xx_contour.jpg';
if n<10
 fname04(6)='0';
 fname04(7)=sn(1);
else
 fname04(6)=sn(1);
 fname04(7)=sn(2);
end

fname05 = 'veli_xx_color.jpg';
if n<10
 fname05(6)='0';
 fname05(7)=sn(1);
else
 fname05(6)=sn(1);
 fname05(7)=sn(2);
end

fname06 = 'veli_xx_contour.jpg';
if n<10
 fname06(6)='0';
 fname06(7)=sn(1);
else
 fname06(6)=sn(1);
 fname06(7)=sn(2);
end

fname07 = 'vele_xx_color.jpg';
if n<10
 fname07(6)='0';
 fname07(7)=sn(1);
else
 fname07(6)=sn(1);
 fname07(7)=sn(2);
end

fname08 = 'vele_xx_contour.jpg';
if n<10
 fname08(6)='0';
 fname08(7)=sn(1);
else
 fname08(6)=sn(1);
 fname08(7)=sn(2);
end

fname09 = 'Bmag_xx_color.jpg';
if n<10
 fname09(6)='0';
 fname09(7)=sn(1);
else
 fname09(6)=sn(1);
 fname09(7)=sn(2);
end

fname10 = 'Bmag_xx_contour.jpg';
if n<10
 fname10(6)='0';
 fname10(7)=sn(1);
else
 fname10(6)=sn(1);
 fname10(7)=sn(2);
end

fname11 = 'Emag_xx_color.jpg';
if n<10
 fname11(6)='0';
 fname11(7)=sn(1);
else
 fname11(6)=sn(1);
 fname11(7)=sn(2);
end

fname12 = 'Emag_xx_contour.jpg';
if n<10
 fname12(6)='0';
 fname12(7)=sn(1);
else
 fname12(6)=sn(1);
 fname12(7)=sn(2);
end


disp([' fname01 = ',fname01]);
disp([' fname02 = ',fname02]);
disp([' fname03 = ',fname03]);
disp([' fname04 = ',fname04]);
disp([' fname05 = ',fname05]);
disp([' fname06 = ',fname06]);
disp([' fname07 = ',fname07]);
disp([' fname08 = ',fname08]);
disp([' fname09 = ',fname09]);
disp([' fname10 = ',fname10]);
disp([' fname11 = ',fname11]);
disp([' fname12 = ',fname12]);

figure(1);  print('-djpeg95',fname01);
figure(2);  print('-djpeg95',fname02);
figure(3);  print('-djpeg95',fname03);
figure(4);  print('-djpeg95',fname04);
figure(5);  print('-djpeg95',fname05);
figure(6);  print('-djpeg95',fname06);
figure(7);  print('-djpeg95',fname07);
figure(8);  print('-djpeg95',fname08);
figure(9);  print('-djpeg95',fname09);
figure(10); print('-djpeg95',fname10);
figure(11); print('-djpeg95',fname11);
figure(12); print('-djpeg95',fname12);