//use std::fs::File;
//use std::io::Write;

fn main(){
    let dz = 1.0;
    let dt = 5.0;
    let ntot = (5.0*24.0*3600.0/dt) as isize;
    let nz: usize = 500;

    let mut u = vec![0.0; nz + 1];
    let mut v = vec![0.0; nz + 1];
    let mut z = vec![0.0; nz];

    let ncfile = "exercise1.nc";

    // Create NetCDF
    create_nc(ncfile, nz);

    {
        // Compute depth levels
        for i in 0..nz {
            let ii: f32 = (i+1) as f32;
            z[i] = -ii*dz+0.5*dz;
        }
        // re-open it in append mode
        let mut nc = netcdf::append(&ncfile).unwrap();
        let var = &mut nc.variable_mut("depth").unwrap();
        let res = var.put_values(&z, None, None);
        assert_eq!(res.unwrap(), ());
        // close it (done when `nc` goes out of scope)
    };


    // Run numerical scheme
    let mut j: usize = 0;
    for i in 1..=ntot {
        let (_u, _v) = dynamics(nz, &mut u, &mut v);
        if i % (24.0*3600.0/dt) as isize == 0 {

            // compute time
            let ii: f32 = i as f32;

            // time
            {
                // re-open it in append mode
                let mut nc = netcdf::append(&ncfile).unwrap();
                let var = &mut nc.variable_mut("time").unwrap();
                let res = var.put_values(&[ii*dt], Some(&[j]), None);
                assert_eq!(res.unwrap(), ());
                // close it (done when `nc` goes out of scope)
            };

            // u eastward
            {
                // re-open it in append mode
                let mut nc = netcdf::append(&ncfile).unwrap();
                let var = &mut nc.variable_mut("u").unwrap();
                let un = u[1..nz+1].to_vec();
                let res = var.put_values(&un, Some(&[j, 0]), None);
                assert_eq!(res.unwrap(), ());
                // close it (done when `nc` goes out of scope)
            };

            j = j + 1;
        };
    }

    /*
    // create txt file
    let mut f = File::create("zuv.dat").expect("unable to create");
    writeln!(&mut f, "{} {} {}", "depth", "u", "v").unwrap();

    for i in 1..=nz{
        let ii: f32 = i as f32;
        let z:f32 = -ii*dz+0.5*dz;
        writeln!(&mut f, "{} {} {}", z.to_string(), u[i].to_string(), v[i].to_string()).unwrap();
    }
    */
}

fn dynamics(nz: usize, u: &mut Vec<f32>, v: &mut Vec<f32>)-> (Vec<f32>, Vec<f32>) {

    let az = 5.0e-2;
    let tx = 0.0;
    let ty = 0.5;
    let rho0 = 1028.0;
    let dz = 1.0;
    let dt = 5.0;

    let alpha = 1.0e-4*dt;
    let beta = 0.25*alpha*alpha;

    let a = alpha;
    let b = 1.0-beta;
    let c = 1.0+beta;

    let mut un = vec![0.0; nz + 1];
    let mut vn = vec![0.0; nz + 1];

    let atop = az; // general: 0.5*(Az[0] + Az[1])
    let abot = az; // general: 0.5*(Az[1] - Az[2])

    // Surface boundary condition
    u[0] = u[1] + dz*tx/rho0/atop;
    v[0] = v[1] + dz*ty/rho0/atop;

    // Numerical scheme
    for i in 1..nz {
        let diffu: f32 = dt*(atop*(u[i-1] - u[i])/dz -abot*(u[i]-u[i+1])/dz)/dz;
        let diffv: f32 = dt*(atop*(v[i-1] - v[i])/dz -abot*(v[i]-v[i+1])/dz)/dz;
        un[i] = (b*u[i] + a*v[i] + 0.5*a*diffv + diffu)/c;
        vn[i] = (b*v[i] - a*u[i] - 0.5*a*diffu + diffv)/c;
    }

    //bottom boundary condition
    un[nz] = un[nz-1];
    vn[nz] = vn[nz-1];

    for i in 1..nz+1{
        u[i] = un[i];
        v[i] = vn[i];
    }

    return (u.to_vec(), v.to_vec())
}

fn create_nc(ncfile: &str, nz: usize){
    // Create a new NetCDF
    let mut nc = netcdf::create(ncfile).unwrap();

    // Dimensions
    nc.add_dimension("depth", nz).unwrap();
    nc.add_unlimited_dimension("time").unwrap();

    // Variables
    let mut depthnc = nc.add_variable::<f32>(
                "depth",
                &["depth"],
    ).unwrap();
    // Metadata can be added to the variable
    depthnc.add_attribute("long_name", "depth").unwrap();
    depthnc.add_attribute("units", "metres").unwrap();
    depthnc.add_attribute("axis", "Z").unwrap();

    let mut tnc = nc.add_variable::<f32>(
                "time",
                &["time"],
    ).unwrap();
    // Metadata can be added to the variable
    tnc.add_attribute("long_name", "time").unwrap();
    tnc.add_attribute("units", "seconds since 2021-06-28 00:00:00").unwrap();
    tnc.add_attribute("calendar", "standard").unwrap();
    tnc.add_attribute("axis", "T").unwrap();

    let mut unc = nc.add_variable::<f32>(
                "u",
                &["time", "depth"],
    ).unwrap();
    // Metadata can be added to the variable
    unc.add_attribute("standard_name", "eastward_sea_water_velocity").unwrap();
    unc.add_attribute("long_name", "Eastward velocity").unwrap();
    unc.add_attribute("units", "m s-1").unwrap();

    let mut vnc = nc.add_variable::<f32>(
                "v",
                &["time", "depth"],
    ).unwrap();
    // Metadata can be added to the variable
    vnc.add_attribute("standard_name", "northward_sea_water_velocity").unwrap();
    vnc.add_attribute("long_name", "Northward velocity").unwrap();
    vnc.add_attribute("units", "m s-1").unwrap();
}
