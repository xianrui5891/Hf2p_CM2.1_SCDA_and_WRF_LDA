begin
    ; Load the WRF input file
    aa = systemfunc("rm -rf output_PS.nc ")
    f_in = addfile("../wrfinput_d01", "r")

    ; Read the original latitude and longitude
    lat_orig = f_in->XLAT_V
    lon_orig = f_in->XLONG_U

    ; 选择新维度
    lat = todouble( lat_orig(0,:,1) )
    lon = todouble( lon_orig(0,1,:) )

    ; 读取温度数据
    PS_orig1 = f_in->PSFC
    PS_orig  = todouble(PS_orig1)

    n_size = dimsizes(PS_orig)
    
    ; 创建新的温度变量
    PS_new = new((/n_size(0), n_size(1)+1, n_size(2)+1/), "double")
    do i = 0, n_size(1)-1
        do j = 0, n_size(2)-1
            PS_new(0,i,j) = PS_orig(0,i,j)
        end do
        ; write loop content
    end do
    do i = 0, n_size(1)-1
            PS_new(0,i,n_size(2)) = PS_orig(0,i,n_size(2)-1)
    end do
    do j = 0, n_size(2)-1
            PS_new(0,n_size(1),j) = PS_orig(0,n_size(1)-1,j)
    end do
    PS_new(0,n_size(1),n_size(2)) = PS_orig(0,n_size(1)-1,n_size(2)-1)

    ; 输出到新的 NetCDF 文件
    fout = addfile("output_PS.nc", "c")
    setfileoption(fout, "DefineMode", True)
    setfileoption("nc", "Format", "NetCDF4Classic")

    ; 定义维度
    dimNames = (/ "z", "lat", "lon" /)
    dimSizes = (/ n_size(0), n_size(1)+1, n_size(2)+1/)
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

    PS_new!0 = "z"
    PS_new!1 = "lat"
    PS_new!2 = "lon"
    ; T&z = 49
    PS_new&lat = lat
    PS_new&lon = lon
    ; PS_new&0 = 49
    ; PS_new&1 = lat_new
    ; PS_ne&2 =  lon_new
        ; 定义变量
    filevardef(fout, "lat", typeof(lat), getvardims(lat))
    filevardef(fout, "lon", typeof(lon), getvardims(lon))
    filevardef(fout, "PS", typeof(PS_new), getvardims(PS_new))

    ; 写入数据
    fout->lat  = lat
    fout->lon  = lon
    fout->time = time
    fout->PS    = PS_new


    setfileoption(fout, "DefineMode", False)
end

