
%-----------------------------------------------------------------------
%dyneye
%Copyright (C) Floris van Breugel, 2013.
%  
%florisvb@gmail.com
%
%Released under the GNU GPL license, Version 3
%
%This file is part of dyneye.
%
%dyneye is free software: you can redistribute it and/or modify it
%under the terms of the GNU General Public License as published
%by the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%    
%dyneye is distributed in the hope that it will be useful, but WITHOUT
%ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
%FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
%License for more details.
%
%You should have received a copy of the GNU General Public
%License along with dyneye.  If not, see <http://www.gnu.org/licenses/>.
%
%------------------------------------------------------------------------


load('data.mat')

xf0 = [10;-1];

f = @(x,u,w) x + ([x(2); u])*dt;
h = @(x,u,w) x(2)/x(1) + w;
[xf, P] = ukf_sqrt(optic_flow', xf0, f, h, Q, R, control');

subplot(1,3,1);
plot(time, position, 'blue', time, xf(1,:), 'red');
subplot(1,3,2);
plot(time, velocity, 'blue', time, xf(2,:), 'red');
subplot(1,3,3);

l = size(P);
Pp = reshape(P(1,1,:),l(3),1);

error = position - xf(1,:)';
plot(time, 3*sqrt(Pp), 'black', time, -3*sqrt(Pp), 'black', time, error, 'red');

savefile = 'ukf_matlab_output.mat'
save(savefile, 'time', 'xf');
