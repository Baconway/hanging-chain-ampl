# ==============================================================================
# Mô hình: xich_chuong_ngai_vat.mod
# Bài tập 4: Bài toán dây xích treo (Câu 3 - Có chướng ngại vật)
# ==============================================================================

# Số lượng mắt xích cứng (m > 1)
param so_mat_xich integer > 1;

# Tổng chiều dài của dây xích (L > 0)
param chieu_dai > 0;

# Tọa độ điểm đầu cuối bên phải: dây xích nối từ gốc (0,0) đến (toa_do_a, toa_do_b)
param toa_do_a;
param toa_do_b;

# Hệ số trọng số thế năng của từng mắt xích: c[k] = m - k + 0.5
param he_so_c {k in 1..so_mat_xich} := so_mat_xich - k + 0.5;

# Các thông số của hàm chướng ngại vật:
# obst(x) = chieu_cao_h * exp( - (x - tam_mu)^2 / phuong_sai ) - do_dich_nu
param chieu_cao_h >= 0;
param tam_mu;
param phuong_sai > 0;
param do_dich_nu;

# Số điểm rời rạc hóa trên mỗi mắt xích để kiểm tra va chạm với chướng ngại vật
param so_diem_chia integer > 0 default 20;

# Biến quyết định: độ tăng tọa độ theo phương ngang và phương đứng của mỗi mắt xích
var x {1..so_mat_xich};
var y {1..so_mat_xich};

# Tọa độ tích lũy của các khớp nối:
# Khớp 0 được cố định tại gốc tọa độ (0, 0)
# Khớp k có tọa độ là (tổng_{i=1}^k x_i, tổng_{i=1}^k y_i)
var toa_do_khop_x {k in 0..so_mat_xich} = if k = 0 then 0.0 else sum {i in 1..k} x[i];
var toa_do_khop_y {k in 0..so_mat_xich} = if k = 0 then 0.0 else sum {i in 1..k} y[i];

# Rời rạc hóa tọa độ dọc theo từng mắt xích k thuộc 1..so_mat_xich:
# Đoạn thẳng nối giữa khớp (k-1) và khớp k được tham số hóa bởi tỷ lệ i / so_diem_chia thuộc [0, 1]
var diem_chia_x {k in 1..so_mat_xich, i in 0..so_diem_chia} =
    (1 - i / so_diem_chia) * toa_do_khop_x[k-1] + (i / so_diem_chia) * toa_do_khop_x[k];
var diem_chia_y {k in 1..so_mat_xich, i in 0..so_diem_chia} =
    (1 - i / so_diem_chia) * toa_do_khop_y[k-1] + (i / so_diem_chia) * toa_do_khop_y[k];

# Hàm mục tiêu: Cực tiểu hóa tổng thế năng trọng trường của dây xích
minimize TheNang:
    sum {k in 1..so_mat_xich} he_so_c[k] * y[k];

# Các ràng buộc vật lý:
# 1. Ràng buộc chiều dài mắt xích cứng bằng nhau: x_k^2 + y_k^2 = (L/m)^2
subject to DoDaiMatXich {k in 1..so_mat_xich}:
    x[k]^2 + y[k]^2 = (chieu_dai / so_mat_xich)^2;

# 2. Ràng buộc tọa độ ngang của điểm cuối bên phải: tổng x_k = a
subject to DiemCuoiX:
    sum {k in 1..so_mat_xich} x[k] = toa_do_a;

# 3. Ràng buộc tọa độ đứng của điểm cuối bên phải: tổng y_k = b
subject to DiemCuoiY:
    sum {k in 1..so_mat_xich} y[k] = toa_do_b;

# 4. Ràng buộc tránh chướng ngại vật: toàn bộ các điểm rời rạc trên dây xích phải nằm trên hoặc tiếp xúc chướng ngại vật
subject to TranhChuongNgaiVat {k in 1..so_mat_xich, i in 0..so_diem_chia}:
    diem_chia_y[k,i] >= chieu_cao_h * exp( - (diem_chia_x[k,i] - tam_mu)^2 / phuong_sai ) - do_dich_nu;

