function spacecraft = theChaserSpacecraft(cx,cy,theta,r,wdock)
    % This function takes as input the position of the center of mass of
    % the chaser spacecraft and its orientation (from the solution), and
    % its width properties, and outputs a list of vector points of its
    % outline. For plotting purposes.
    
    % Define a rotation matrix based on the chaser orientation
    rotmat = [cos(theta), -sin(theta); 
              sin(theta), cos(theta)];
          
    % Define the points that make-up the outline of the chaser spacecraft
    sc = [cx,   cy;
          cx+r, cy;
          cx+r, cy+r;
          cx-r, cy+r;
          cx-r, cy-r;
          cx+r, cy-r;
          cx+r, cy-r/3;
          cx+r+wdock, cy;
          cx+r, cy+r/3;
          cx+r, cy-r/3;
          cx+r, cy];

    % Transform the spacecraft to be at the origin in order to apply the
    % rotation matrix
    sccenter = [ones(1,length(sc))*cx; ones(1,length(sc))*cy].' ;
    scbf = sc - sccenter;
    
    % Apply the rotation matrix to each of the points
    scbf_rot = zeros(size(scbf));
    for i=1:length(sc)
        scbf_rot(i,:) = rotmat*scbf(i,:).';
    end
    
    % Transform the spacecraft back to its correct position at cx,cy
    spacecraft = scbf_rot + sccenter;
end