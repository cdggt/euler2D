function pt = project_state_to_sphere(domain,origin, pt, rad)
    
pt = origin+rad*util.hat(domain,pt-origin);

end