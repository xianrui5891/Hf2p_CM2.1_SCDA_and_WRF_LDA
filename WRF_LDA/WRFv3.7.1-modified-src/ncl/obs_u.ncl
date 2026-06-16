begin
    ; Load the WRF input file
    aa = systemfunc("rm -rf output_U.nc ")
    f_in = addfile("../wrfinput_d01", "r")

    ; Read the original latitude and longitude
    lat_orig = f_in->XLAT_V
    lon_orig = f_in->XLONG_U

    ; 选择新维度
    lat = todouble( lat_orig(0,:,1) )
    lon = todouble( lon_orig(0,1,:) )

    ; 读取温度数据
    U_orig1 = f_in->U(0,:,:,:)
    U_orig  = todouble(U_orig1)

    n_size = dimsizes(U_orig)
    
    ; 创建新的温度变量
    U_new = new((/n_size(0), n_size(1)+1, n_size(2)/), "double")
    do i = 0, n_size(1)-1
        do j = 0, n_size(2)-1
            do k = 0, n_size(0)-1
                U_new(k,i,j) = U_orig(k,i,j)
            end do
        end do
    end do
    ; do i = 0, n_size(1)-1
    ;     do k = 0, n_size(0)-1
    ;         U_new(k,i,n_size(2)) = U_orig(k,i,n_size(2)-1)
    ;     end do
    ; end do
    do j = 0, n_size(2)-1
        do k = 0, n_size(0)-1
            U_new(k,n_size(1),j) = U_orig(k,n_size(1)-1,j)
        end do
    end do
    U_new(n_size(0)-1,n_size(1),n_size(2)-1) = U_orig(n_size(0)-1,n_size(1)-1,n_size(2)-1)



    ; 输出到新的 NetCDF 文件
    fout = addfile("output_U.nc", "c")
    setfileoption(fout, "DefineMode", True)
    setfileoption("nc", "Format", "NetCDF4Classic")

    ; 定义维度
    dimNames = (/ "z", "lat", "lon" /)
    dimSizes = (/ n_size(0), n_size(1)+1, n_size(2) /)
    dimUnlim = (/ False, False, False /)
    filedimdef(fout, dimNames, dimSizes, dimUnlim)
    time = todouble(1.0) 
  


    ; 设置变量的属性
    lat!0 = "lat"
    lat&lat = lat 
    lon!0 = "lon"
    lon&lon = lon
    time!0 = "time"
    time&time = time

    U_new!0 = "z"
    U_new!1 = "lat"
    U_new!2 = "lon"
    ; T&z = 49
    U_new&lat = lat
    U_new&lon = lon
    ; U_new&0 = 49
    ; U_new&1 = lat_new
    ; U_ne&2 =  lon_new
        ; 定义变量
    filevardef(fout, "lat", typeof(lat), getvardims(lat))
    filevardef(fout, "lon", typeof(lon), getvardims(lon))
    filevardef(fout, "U", typeof(U_new), getvardims(U_new))

    ; 写入数据
    fout->lat  = lat
    fout->lon  = lon
    fout->time = time
    fout->U    = U_new


    setfileoption(fout, "DefineMode", False)
end

