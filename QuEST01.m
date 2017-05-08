function [q_BO_est] = QuEST01(v_mag_earth,v_grav_earth,v_mag_mobile,v_grav_mobile, weight)
    % q_BO_est = Estimated quaternion
    % b1,b2 = the 2 measurement vectors in body frames
    % r1,r2 = the 2 measurement vectors in reference frames (orbit/inertial; orbit used)
    % b1 is the measurement with very high accurace as compared to b2
    % a1 and a2 are Relative weights of b1 and b2

    %%%% start TRIAD %%%%
    b1 = v_mag_earth;
    b2 = v_grav_earth;
    r1 = v_mag_mobile;
    r2 = v_grav_mobile;

    a1 = weight;
    a2 = 1-a1;

    %b3
    b3 = cross(b1,b2);
    b3 = b3/norm(b3) ;
    %r3
    r3 = cross(r1,r2);
    r3 = r3/norm(r3) ;

    %mu
    mu = ( 1 + dot(b3,r3) ) * ( a1*dot(b1,r1) + a2*dot(b2,r2) ) ;
    mu =  mu + dot( cross(b3,r3) , ( a1*cross(b1,r1) + a2*cross(b2,r2) ) ) ;

    %nu
    nu = dot( ( b3 + r3 ) , ( a1*cross(b1,r1) + a2*cross(b2,r2) )  ) ;

    %rho
    rho = sqrt( mu*mu + nu*nu );

    %q
    if( mu > 0)
        k_q_triad = (1/( 2*sqrt( rho*(rho+mu)*( 1 + dot(b3,r3) ) ) ) );
        v_q_triad = ( rho + mu )*( cross(b3,r3) ) + nu*(b3 + r3) ;
        s_q_triad = ( rho + mu )* ( 1 + dot(b3,r3) ) ;
        q_triad = [v_q_triad ; s_q_triad];
        q_triad = k_q_triad * q_triad ;
    else
        k_q_triad = (1/( 2*sqrt( rho*(rho-mu)*( 1 + dot(b3,r3) ) ) ) );
        v_q_triad = ( nu )*( cross(b3,r3) ) + ( rho - mu )*(b3 + r3) ;
        s_q_triad = ( nu )* ( 1 + dot(b3,r3) ) ;
        q_triad = [v_q_triad ; s_q_triad];
        q_triad = k_q_triad * q_triad ;    
    end

    q_triad = q_triad / sqrt(dot(q_triad,q_triad)) ;
    q_BO_est = q_triad;

