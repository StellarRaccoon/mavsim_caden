classdef mav_viewer < handle
    %
    %    Create spacecraft animation
    %
    %--------------------------------
    properties
        body_handle
    	Vertices
    	Faces
    	facecolors
        plot_initialized
    end
    %--------------------------------
    methods
        %------constructor-----------
        function self = mav_viewer
            self.body_handle = [];
            [self.Vertices, self.Faces, self.facecolors] = self.define_mav();
            self.plot_initialized = 0;           
        end
        %---------------------------
        function self=update(self, state)
            if self.plot_initialized==0
                figure(1); clf;
                self=self.drawBody(state.pn, state.pe, -state.h,...
                                    state.phi, state.theta, state.psi);
                title('MAV')
                xlabel('East')
                ylabel('North')
                zlabel('-Down')
                view(32,47)  % set the vieew angle for figure
                % axis([-10,10,-10,10,-10,10]);
                axis([-30,30,0,200,70,120]);
                axis equal;
                hold on
                grid on
                self.plot_initialized = 1;
            else
                self=self.drawBody(state.pn, state.pe, -state.h,... 
                                state.phi, state.theta, state.psi);

            end
        end
        %---------------------------
        function self = drawBody(self, pn, pe, pd, phi, theta, psi)
            Vertices = self.rotate(self.Vertices, phi, theta, psi);   % rotate rigid body  
            Vertices = self.translate(Vertices, pn, pe, pd);     % translate after rotation
            % transform vertices from NED to ENU (for matlab rendering)
            R = [...
                0, 1, 0;...
                1, 0, 0;...
                0, 0, -1;...
                ];
            Vertices = R*Vertices;
            if isempty(self.body_handle)
                self.body_handle = patch('Vertices', Vertices', 'Faces', self.Faces,...
                                        'FaceVertexCData',self.facecolors,...
                                        'FaceColor','flat');
            else
                set(self.body_handle,'Vertices',Vertices','Faces',self.Faces);
                drawnow
            end
        end 
        %---------------------------
        function pts=rotate(self, pts, phi, theta, psi)
            % define rotation matrix (right handed)
            R_roll = [...
                        1, 0, 0;
                        0, cos(phi), sin(phi);
                        0, -sin(phi), cos(phi)];
            R_pitch = [...
                        cos(theta), 0, -sin(theta);
                        0, 1, 0;
                        sin(theta), 0, cos(theta)];
            R_yaw = [...
                        cos(psi), sin(psi), 0;
                        -sin(psi), cos(psi), 0;
                        0, 0, 1];
            R = R_roll*R_pitch*R_yaw;   % inertial to body
            R = R';  % body to inertial
            % rotate vertices
            pts = R*pts;
        end
        %---------------------------
        % translate vertices by pn, pe, pd
        function pts = translate(self, pts, pn, pe, pd)
            pts = pts + repmat([pn;pe;pd],1,size(pts,2));
        end
        %---------------------------
        function [V, F, colors] = define_mav(self)
            fuse_h = 3;
            fuse_w = 2;
            fuse_l1 = 5;
            fuse_l2 = 2;
            fuse_l3 = 10;
            wing_l = 3;
            wing_w = 10;
            tail_h = 3;
            tailwing_l = 2;
            tailwing_w = 3;
            % Define the vertices (physical location of vertices)
            V = [
                fuse_l1, 0, 0;                              % point 1
                fuse_l2, fuse_w/2, -fuse_h/2;               % point 2
                fuse_l2, -fuse_w/2, -fuse_h/2;              % point 3
                fuse_l2, -fuse_w/2, fuse_h/2;               % point 4
                fuse_l2, fuse_w/2, fuse_h/2;                % point 5
                -fuse_l3, 0, 0;                             % point 6
                0, wing_w/2, 0;                             % point 7
                -wing_l, wing_w/2, 0;                       % point 8
                -wing_l, -wing_w/2, 0;                      % point
                0, -wing_w/2, 0;                            % point 10
                -fuse_l3+tailwing_l, tailwing_w/2, 0;       % point 11
                -fuse_l3, tailwing_w/2, 0;                  % point 12
                -fuse_l3, -tailwing_w/2, 0;                 % point 13
                -fuse_l3+tailwing_l, -tailwing_w/2, 0;      % point 14
                -fuse_l3+tailwing_l, 0, 0;                  % point 15
                -fuse_l3, 0, -tail_h;                       % point 16
            ]';

            % define faces as a list of vertices numbered above
            F = [...
                    1, 2, 3;        % top nose
                    1, 3, 4;        % port nose
                    1, 4, 5;        % bottom nose 
                    1, 5, 2;        % starboard nose
                    3, 6, 4;        % port fusalaage
                    2, 3, 6;        % top fusalage
                    2, 5, 6;        % starboard fusalage
                    4, 5, 6;        % bottom fusalage 
                    7, 8, 9;        % wing 1
                    7, 9, 10;       % wing 2
                    11, 12, 13;     % Horiz stabalizer 1
                    11, 13, 14;     % h stabalizer 2
                    6, 15, 16;      % v stabalizer
                ];

            % define colors for each face    
            myred = [1, 0, 0];
            mygreen = [0, 1, 0];
            myblue = [0, 0, 1];
            myyellow = [1, 1, 0];
            mycyan = [0, 1, 1];

            colors = [...
                myred;      % top nose
                mygreen;    % port nose
                myblue;     % bottom nose 
                myyellow;   % starboard nose
                myred;      % port fusalaage
                mygreen;    % top fusalage
                myblue;     % starboard fusalage
                myyellow;   % bottom fusalage 
                mycyan;     % wing 1
                myyellow;   % wing 2
                mycyan;     % Horiz stabalizer 1
                myyellow;   % h stabalizer 2
                myyellow;   % v stabalizer
            ];
        end
    end
end