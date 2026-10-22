# ==============================================================================
# Mô hình: xich_co_ban.mod
# Bài tập 4: Bài toán dây xích treo (Câu 2 - Mô hình cơ bản với 3 mắt xích)
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

# Biến quyết định: độ tăng tọa độ theo phương ngang và phương đứng của mỗi mắt xích
var x {1..so_mat_xich};
var y {1..so_mat_xich};

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

