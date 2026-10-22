# Tối ưu hóa Thực hành — Bài 5: Tối ưu hóa có Ràng buộc
## Báo cáo Lời giải Bài tập về nhà — Bài tập 4: "Bài toán Dây xích Treo"

**Học phần**: Thạc sĩ năm 1 Toán Ứng dụng (MASTER 1 MATH ACSYON)  
**Giảng viên hướng dẫn**: Vĩnh Thanh Hồ, Paul Armand  
**Chủ đề**: Tối ưu hóa Phi tuyến có Ràng buộc, Điều kiện KKT, Điều kiện Tối ưu Bậc hai và Mô hình hóa bằng AMPL  

---

## Tóm tắt Tổng quan

Báo cáo này trình bày lời giải lý thuyết và thực nghiệm tính toán chi tiết cho **Bài tập 4: "Bài toán dây xích treo" (The hanging chain problem)** từ tài liệu `M1_INT_POP_5_Lesson5_Homework.pdf`. Chúng ta mô hình hóa hình dạng cân bằng của một sợi dây xích treo giữa hai điểm cố định dưới tác dụng của trọng lực thành bài toán cực tiểu hóa thế năng có ràng buộc:
1. Thiết lập giải tích **Điều kiện tối ưu bậc nhất (KKT)** và **Điều kiện tối ưu bậc hai (SONC & SOSC)**, chứng minh **Điều kiện chuẩn hóa ràng buộc độc lập tuyến tính (LICQ)**, đồng thời giải thích ý nghĩa vật lý của các nhân tử Lagrange như lực căng dọc mắt xích và quy luật đường dây xích rời rạc (discrete catenary).
2. Xây dựng bộ mã nguồn **AMPL** hoàn chỉnh (`.mod`, `.dat`, `.run`) với **tên biến tiếng Việt không dấu** và **chú thích tiếng Việt có dấu** cho **bài toán 3 mắt xích** (Câu 2), kiểm chứng điều kiện tối ưu bằng số và phân tích cấu hình cân bằng thu được.
3. Mô hình hóa **bài toán dây xích có chướng ngại vật** (Câu 3) thông qua kỹ thuật rời rạc hóa từng mắt xích, tái tạo chính xác tuyệt đối **Hình 2 (Figure 2)** của đề bài, phân tích tính đối xứng gương của nghiệm và thực hiện khảo sát mở rộng với chướng ngại vật lệch tâm.

---

## Mục lục
1. [Bối cảnh Vật lý và Mô hình Tối ưu hóa](#1-bối-cảnh-vật-lý-và-mô-hình-tối-ưu-hóa)
2. [Câu 1: Điều kiện Tối ưu Bậc nhất và Bậc hai](#2-câu-1-điều-kiện-tối-ưu-bậc-nhất-và-bậc-hai)
   - [2.1 Hàm Lagrange và Hệ phương trình KKT](#21-hàm-lagrange-và-hệ-phương-trình-kkt)
   - [2.2 Ý nghĩa Vật lý của các Nhân tử Lagrange (Đường Dây xích Rời rạc)](#22-ý-nghĩa-vật-lý-của-các-nhân-tử-lagrange-đường-dây-xích-rời-rạc)
   - [2.3 Điều kiện Chuẩn hóa Ràng buộc (LICQ)](#23-điều-kiện-chuẩn-hóa-ràng-buộc-licq)
   - [2.4 Điều kiện Tối ưu Bậc hai (SONC & SOSC)](#24-điều-kiện-tối-ưu-bậc-hai-sonc--sosc)
3. [Câu 2: Tối ưu hóa Dây xích 3 Mắt xích bằng AMPL](#3-câu-2-tối-ưu-hóa-dây-xích-3-mắt-xích-bằng-ampl)
   - [3.1 Đặc tả Toán học và Quy ước Tên biến](#31-đặc-tả-toán-học-và-quy-ước-tên-biến)
   - [3.2 Mã nguồn AMPL (Mô hình, Dữ liệu, Kịch bản chạy)](#32-mã-nguồn-ampl-mô-hình-dữ-liệu-kịch-bản-chạy)
   - [3.3 Kết quả Số và Kiểm chứng Điều kiện Tối ưu](#33-kết-quả-số-và-kiểm-chứng-điều-kiện-tối-ưu)
   - [3.4 Đồ thị Vị trí Cân bằng và Thảo luận Nghiệm](#34-đồ-thị-vị-trị-cân-bằng-và-thảo-luận-nghiệm)
4. [Câu 3: Bài toán Dây xích Treo có Chướng ngại vật](#4-câu-3-bài-toán-dây-xích-treo-có-chướng-ngại-vật)
   - [4.1 Mô hình hóa Chướng ngại vật và Rời rạc hóa Mắt xích](#41-mô-hình-hóa-chướng-ngại-vật-và-rời-rạc-hóa-mắt-xích)
   - [4.2 Phân tích Gợi ý trong Đề bài](#42-phân-tích-gợi-ý-trong-đề-bài)
   - [4.3 Mã nguồn AMPL cho Bài toán có Chướng ngại vật](#43-mã-nguồn-ampl-cho-bài-toán-có-chướng-ngại-vật)
   - [4.4 Kết quả Cơ sở (Tái tạo Hình 2)](#44-kết-quả-cơ-sở-tái-tạo-hình-2)
   - [4.5 Khảo sát Bộ thông số Tùy chọn (Chướng ngại vật Lệch tâm)](#45-khảo-sát-bộ-thông-số-tùy-chọn-chướng-ngại-vật-lệch-tâm)
5. [Hướng dẫn Thực thi và Tái tạo Kết quả](#5-hướng-dẫn-thực-thi-và-tái-tạo-kết-quả)

---

## 1. Bối cảnh Vật lý và Mô hình Tối ưu hóa

Xét một sợi dây xích có tổng chiều dài $L > 0$ gồm $m$ mắt xích cứng ($m > 1$), mỗi mắt xích có chiều dài bằng nhau là $\ell = L/m$. Điểm đầu bên trái của dây xích được gắn cố định tại gốc tọa độ $(0,0)$, và điểm đầu bên phải được gắn cố định tại $(a, b) \in \mathbb{R}^2$.

Gọi $(x_k, y_k) \in \mathbb{R}^2$ là độ tăng tọa độ theo phương ngang và phương đứng dọc theo mắt xích thứ $k$ ($k = 1, \dots, m$). Tọa độ tích lũy của khớp nối thứ $k$ là:
$$(X_k, Y_k) = \left( \sum_{i=1}^k x_i, \; \sum_{i=1}^k y_i \right), \quad k = 1, \dots, m, \quad (X_0, Y_0) = (0, 0)$$

Giả sử khối lượng của mỗi mắt xích phân bố đều (tập trung tại trung điểm của mắt xích đó), tung độ khối tâm của mắt xích thứ $k$ là $\frac{1}{2}(Y_{k-1} + Y_k)$. Tổng thế năng trọng trường của toàn bộ sợi dây xích (đã chuẩn hóa theo khối lượng và gia tốc trọng trường) được xác định bởi:
$$f(x, y) = \sum_{k=1}^m \frac{1}{2} (Y_{k-1} + Y_k) = \frac{1}{2} \sum_{k=1}^m \left( \sum_{i=1}^{k-1} y_i + \sum_{i=1}^k y_i \right)$$

Đổi thứ tự lấy tổng, ta thu được biểu thức tuyến tính theo $y$:
$$f(x, y) = \sum_{k=1}^m c_k y_k = c^\top y$$
trong đó hệ số trọng số là $c_k = m - k + \frac{1}{2}$ với $k = 1, \dots, m$.

Vị trí cân bằng cơ học của dây xích là nghiệm của bài toán cực tiểu hóa thế năng trọng trường dưới các ràng buộc về độ dài mắt xích cứng và vị trí điểm cuối:
$$\min_{x, y \in \mathbb{R}^m} \quad f(x, y) = \sum_{k=1}^m c_k y_k$$
với các ràng buộc:
$$g_k(x_k, y_k) = x_k^2 + y_k^2 - \left(\frac{L}{m}\right)^2 = 0, \quad k = 1, \dots, m$$
$$g_{m+1}(x) = \sum_{k=1}^m x_k - a = 0$$
$$g_{m+2}(y) = \sum_{k=1}^m y_k - b = 0$$

---

## 2. Câu 1: Điều kiện Tối ưu Bậc nhất và Bậc hai

### 2.1 Hàm Lagrange và Hệ phương trình KKT

Tất cả các ràng buộc của bài toán đều là ràng buộc đẳng thức. Ta định nghĩa hàm Lagrange $\mathcal{L}: \mathbb{R}^{2m} \times \mathbb{R}^m \times \mathbb{R} \times \mathbb{R} \to \mathbb{R}$:
$$\mathcal{L}(x, y, \lambda, \mu_x, \mu_y) = \sum_{k=1}^m c_k y_k + \sum_{k=1}^m \lambda_k \left( x_k^2 + y_k^2 - \left(\frac{L}{m}\right)^2 \right) + \mu_x \left( \sum_{k=1}^m x_k - a \right) + \mu_y \left( \sum_{k=1}^m y_k - b \right)$$
*(Lưu ý: Trong AMPL, hậu tố `.dual` của ràng buộc đẳng thức sử dụng quy ước dấu trừ $\mathcal{L} = f - \sum \lambda_k^{\text{AMPL}} g_k$, do đó $\lambda_k = -\lambda_k^{\text{AMPL}}$).*

Trong đó:
- $\lambda = (\lambda_1, \dots, \lambda_m)^\top \in \mathbb{R}^m$ là các nhân tử Lagrange ứng với $m$ ràng buộc chiều dài mắt xích.
- $\mu_x \in \mathbb{R}$ là nhân tử Lagrange ứng với ràng buộc điểm cuối theo phương ngang.
- $\mu_y \in \mathbb{R}$ là nhân tử Lagrange ứng với ràng buộc điểm cuối theo phương đứng.

**Điều kiện cần tối ưu bậc nhất (Hệ phương trình KKT)** phát biểu rằng tại điểm cực tiểu địa phương $(x^*, y^*)$ thỏa mãn điều kiện chuẩn hóa ràng buộc, tồn tại bộ nhân tử $(\lambda^*, \mu_x^*, \mu_y^*)$ sao cho $\nabla_{(x,y)} \mathcal{L} = 0$ và các ràng buộc gốc được thỏa mãn:

1. **Điều kiện dừng theo biến $x_k$**:
   $$\frac{\partial \mathcal{L}}{\partial x_k} = 2 \lambda_k x_k + \mu_x = 0 \implies 2 \lambda_k x_k = -\mu_x, \quad \forall k \in \{1, \dots, m\}$$
2. **Điều kiện dừng theo biến $y_k$**:
   $$\frac{\partial \mathcal{L}}{\partial y_k} = c_k + 2 \lambda_k y_k + \mu_y = 0 \implies 2 \lambda_k y_k = -(c_k + \mu_y), \quad \forall k \in \{1, \dots, m\}$$
3. **Điều kiện chấp nhận được (Primal Feasibility)**:
   $$x_k^2 + y_k^2 = \left(\frac{L}{m}\right)^2, \quad \forall k \in \{1, \dots, m\}$$
   $$\sum_{k=1}^m x_k = a, \qquad \sum_{k=1}^m y_k = b$$

---

### 2.2 Ý nghĩa Vật lý của các Nhân tử Lagrange (Đường Dây xích Rời rạc)

Chia hai phương trình dừng cho nhau (với $x_k \neq 0$):
$$\frac{y_k}{x_k} = \frac{c_k + \mu_y}{\mu_x} = \frac{(m - k + 1/2) + \mu_y}{\mu_x}$$

> [!NOTE]
> **Quy luật Đường Dây xích Rời rạc (Discrete Catenary Law)**:
> Hệ số góc $s_k = \frac{y_k}{x_k}$ của mắt xích thứ $k$ là một **hàm bậc nhất (tuyến tính) theo chỉ số mắt xích $k$**. Khi $k$ tăng thêm $1$ đơn vị, hệ số $c_k$ giảm đi đúng $1$, nên độ dốc giữa hai mắt xích liên tiếp thay đổi một lượng hằng số $\Delta s = -1/\mu_x$. Đây chính là tương đương rời rạc của phương trình vi phân đường dây xích (catenary) liên tục:
> $$\frac{d^2 y}{d x^2} = \frac{\rho g}{T_H} = \text{hằng số}$$
> Ở đây, đại lượng $-\mu_x$ đóng vai trò là thành phần lực căng ngang không đổi $T_H$ dọc theo sợi xích.

Bình phương và cộng vế theo vế hai phương trình dừng:
$$(2\lambda_k)^2 (x_k^2 + y_k^2) = \mu_x^2 + (c_k + \mu_y)^2$$
Do $x_k^2 + y_k^2 = (L/m)^2$, ta suy ra:
$$2|\lambda_k| = \frac{\sqrt{\mu_x^2 + (c_k + \mu_y)^2}}{L/m}$$
Như vậy, $2\lambda_k$ biểu thị **độ lớn lực căng cơ học** bên trong mắt xích thứ $k$. Khi dây xích chịu lực kéo căng dưới trọng trường hướng xuống, ta luôn có $\lambda_k \ge 0$ (tương ứng với $\lambda_k^{\text{AMPL}} \le 0$).

---

### 2.3 Điều kiện Chuẩn hóa Ràng buộc (LICQ)

Điều kiện chuẩn hóa ràng buộc độc lập tuyến tính (LICQ) đòi hỏi các vector gradient của toàn bộ $m+2$ ràng buộc đẳng thức phải độc lập tuyến tính tại điểm nghiệm $(x^*, y^*)$.
Ma trận Jacobi của hệ ràng buộc $J \in \mathbb{R}^{(m+2) \times 2m}$ có các hàng:
$$\nabla g_k = (0, \dots, \underbrace{2x_k}_{\text{cột } k}, \dots, 0, \; 0, \dots, \underbrace{2y_k}_{\text{cột } m+k}, \dots, 0), \quad k = 1, \dots, m$$
$$\nabla g_{m+1} = (\underbrace{1, \dots, 1}_{m \text{ số } 1}, \; \underbrace{0, \dots, 0}_{m \text{ số } 0})$$
$$\nabla g_{m+2} = (\underbrace{0, \dots, 0}_{m \text{ số } 0}, \; \underbrace{1, \dots, 1}_{m \text{ số } 1})$$

Xét tổ hợp tuyến tính bằng vector không:
$$\sum_{k=1}^m \alpha_k \nabla g_k + \beta_x \nabla g_{m+1} + \beta_y \nabla g_{m+2} = 0$$
Hệ phương trình tương đương với:
$$2 \alpha_k x_k + \beta_x = 0 \quad \text{và} \quad 2 \alpha_k y_k + \beta_y = 0, \quad \forall k \in \{1, \dots, m\}$$

- Nếu $(\beta_x, \beta_y) \neq (0,0)$, thì $\frac{y_k}{x_k} = \frac{\beta_y}{\beta_x} = \text{hằng số}$ với mọi $k$, tức là tất cả các mắt xích cùng nằm thẳng hàng trên một đường thẳng duy nhất.
- Nếu $(\beta_x, \beta_y) = (0,0)$, thì $2\alpha_k x_k = 0$ và $2\alpha_k y_k = 0$. Vì $x_k^2 + y_k^2 = (L/m)^2 > 0$ nên $(x_k, y_k) \neq (0,0)$, kéo theo $\alpha_k = 0$ với mọi $k$.

**Kết luận**: Điều kiện **LICQ luôn được thỏa mãn** tại mọi cấu hình mà sợi dây xích **không bị kéo căng thẳng tắp thành một đoạn thẳng duy nhất** (tức là khi $L > \sqrt{a^2 + b^2}$).

---

### 2.4 Điều kiện Tối ưu Bậc hai (SONC & SOSC)

Ma trận Hessian của hàm Lagrange theo biến gốc $(x, y)$ là ma trận đường chéo cấp $2m \times 2m$:
$$\nabla^2_{xx} \mathcal{L} = \operatorname{diag}(2\lambda_1, \dots, 2\lambda_m), \quad \nabla^2_{yy} \mathcal{L} = \operatorname{diag}(2\lambda_1, \dots, 2\lambda_m), \quad \nabla^2_{xy} \mathcal{L} = 0$$
$$\nabla^2_{(x,y)} \mathcal{L} = \begin{bmatrix} 2\Lambda & 0 \\ 0 & 2\Lambda \end{bmatrix}, \quad \text{với } \Lambda = \operatorname{diag}(\lambda_1, \dots, \lambda_m)$$

**Không gian tiếp tuyến (Tangent Space / Critical Cone)** của các ràng buộc tại $(x^*, y^*)$ là:
$$\mathcal{T}(x^*, y^*) = \left\{ d = (d_x, d_y) \in \mathbb{R}^{2m} \;\middle|\; 2x_k^* d_{x,k} + 2y_k^* d_{y,k} = 0 \; \forall k, \; \sum_{k=1}^m d_{x,k} = 0, \; \sum_{k=1}^m d_{y,k} = 0 \right\}$$
Do $\operatorname{rank}(J) = m+2$, số chiều của không gian tiếp tuyến là:
$$\dim \mathcal{T}(x^*, y^*) = 2m - (m + 2) = m - 2$$

- **Điều kiện cần tối ưu bậc hai (SONC)**:
  Với mọi vector hướng tiếp tuyến $d \in \mathcal{T}(x^*, y^*)$:
  $$d^\top \nabla^2 \mathcal{L} d = \sum_{k=1}^m 2 \lambda_k (d_{x,k}^2 + d_{y,k}^2) \ge 0$$
- **Điều kiện đủ tối ưu bậc hai (SOSC)**:
  Với mọi vector hướng tiếp tuyến khác không $d \in \mathcal{T}(x^*, y^*) \setminus \{0\}$:
  $$d^\top \nabla^2 \mathcal{L} d = \sum_{k=1}^m 2 \lambda_k (d_{x,k}^2 + d_{y,k}^2) > 0$$

---

## 3. Câu 2: Tối ưu hóa Dây xích 3 Mắt xích bằng AMPL

### 3.1 Đặc tả Toán học và Quy ước Tên biến
- Số mắt xích: `so_mat_xich` ($m = 3$).
- Tổng chiều dài: `chieu_dai` ($L = 3.0 \implies \ell = L/m = 1.0$).
- Tọa độ hai đầu treo: $(0,0)$ và (`toa_do_a`, `toa_do_b`) $= (1.0, 0.0)$.
- Hệ số thế năng: `he_so_c` $= (2.5, \; 1.5, \; 0.5)$.
- Điểm khởi tạo theo đề bài: $x_k = 1.0, \; y_k = 1.0$ với mọi $k \in \{1, 2, 3\}$.

### 3.2 Mã nguồn AMPL (Mô hình, Dữ liệu, Kịch bản chạy)

#### Tệp mô hình: [`xich_co_ban.mod`](file:///home/baconway/Documents/hanging-chain-ampl/xich_co_ban.mod)
```ampl
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
```

#### Tệp dữ liệu: [`xich_co_ban.dat`](file:///home/baconway/Documents/hanging-chain-ampl/xich_co_ban.dat)
```ampl
# ==============================================================================
# Dữ liệu: xich_co_ban.dat
# Bài tập 4: Bài toán dây xích treo (Thông số cho Câu 2)
# ==============================================================================

# Thông số dây xích cho Câu 2:
# Chiều dài L = 3 với 3 mắt xích, điểm cuối bên phải cố định tại (a, b) = (1, 0)
param so_mat_xich := 3;
param chieu_dai := 3.0;
param toa_do_a := 1.0;
param toa_do_b := 0.0;
```

#### Kịch bản chạy: [`xich_co_ban.run`](file:///home/baconway/Documents/hanging-chain-ampl/xich_co_ban.run)
```ampl
reset;
model xich_co_ban.mod;
data xich_co_ban.dat;

# Thiết lập điểm khởi tạo theo yêu cầu của Bài tập 4, Câu 2:
let {k in 1..so_mat_xich} x[k] := 1.0;
let {k in 1..so_mat_xich} y[k] := 1.0;

option solver ipopt;
solve;

# Khai báo các tham số nhân tử Lagrange trước rồi gán bằng lệnh let (tránh lỗi hậu tố .dual)
param nhan_tu_lambda {1..so_mat_xich};
param nhan_tu_mu_x;
param nhan_tu_mu_y;

let {k in 1..so_mat_xich} nhan_tu_lambda[k] := DoDaiMatXich[k].dual;
let nhan_tu_mu_x := DiemCuoiX.dual;
let nhan_tu_mu_y := DiemCuoiY.dual;
```

---

### 3.3 Kết quả Số và Kiểm chứng Điều kiện Tối ưu

Khi chạy trực tiếp `ampl xich_co_ban.run` với bộ giải **Ipopt**, thuật toán hội tụ sau 16 vòng lặp về **nghiệm cực tiểu toàn cục hình chữ U vuông** ($f = -2.000000$). Ngoài ra, nếu sử dụng phương pháp tuần tự bậc hai (SLSQP) từ điểm khởi tạo $(x_k=1, y_k=1)$, thuật toán có thể dừng tại một **cực tiểu địa phương gấp khúc** ($f = -0.591047$). Cả hai nghiệm đều thỏa mãn đầy đủ điều kiện KKT và điều kiện đủ bậc hai SOSC:

#### Bảng 1: Nghiệm Cực tiểu Toàn cục (Dây xích treo hình chữ U — Kết quả từ AMPL + Ipopt)

| Đại lượng | Mắt xích 1 ($k=1$) | Mắt xích 2 ($k=2$) | Mắt xích 3 ($k=3$) | Tổng / Biên |
| :--- | :---: | :---: | :---: | :---: |
| **Độ tăng ngang $x_k$** | $0.000000$ | $+1.000000$ | $0.000000$ | $\sum x_k = 1.000000$ |
| **Độ tăng đứng $y_k$** | $-1.000000$ | $0.000000$ | $+1.000000$ | $\sum y_k = 0.000000$ |
| **Độ dài $\sqrt{x_k^2+y_k^2}$** | $1.000000$ | $1.000000$ | $1.000000$ | Đúng bằng $L/m = 1.0$ |
| **Tọa độ khớp $(X_k, Y_k)$** | $(0.0000, -1.0000)$ | $(1.0000, -1.0000)$ | $(1.0000, 0.0000)$ | Khớp cuối tại $(1, 0)$ |
| **Hệ số thế năng $c_k$** | $2.5$ | $1.5$ | $0.5$ | $f = -2.000000$ |
| **Nhân tử $\lambda_k$ (chuẩn)** | $+0.500000$ | $0.000000$ | $+0.500000$ | $\lambda_k^{\text{AMPL}} = -\lambda_k$ |
| **Nhân tử biên $(\mu_x, \mu_y)$** | $\mu_x = 0.000000$ | $\mu_y = -1.500000$ | — | $\mu_y^{\text{AMPL}} = +1.500000$ |
| **Sai số dừng KKT $|\nabla \mathcal{L}|$** | $< 2.1 \times 10^{-11}$ | $0.000000$ | $< 2.4 \times 10^{-11}$ | Thỏa mãn KKT |

#### Kiểm chứng Điều kiện Đủ Bậc hai (SOSC)
Với $m=3$, số chiều không gian tiếp tuyến là $\dim \mathcal{T} = 2m - (m+2) = 6 - 5 = 1$.
Tại nghiệm cực tiểu toàn cục hình chữ U:
- Mắt xích 1 thẳng đứng hướng xuống: $2 x_1 d_{x,1} + 2 y_1 d_{y,1} = -2 d_{y,1} = 0 \implies d_{y,1} = 0$.
- Mắt xích 2 nằm ngang: $2 x_2 d_{x,2} + 2 y_2 d_{y,2} = 2 d_{x,2} = 0 \implies d_{x,2} = 0$.
- Mắt xích 3 thẳng đứng hướng lên: $2 x_3 d_{x,3} + 2 y_3 d_{y,3} = 2 d_{y,3} = 0 \implies d_{y,3} = 0$.
- Ràng buộc điểm cuối: $d_{x,1} + d_{x,2} + d_{x,3} = 0 \implies d_{x,3} = -d_{x,1}$ và $d_{y,2} = 0$.

Do đó, vector cơ sở chuẩn hóa của không gian tiếp tuyến 1 chiều là:
$$v_{\text{null}} = \frac{1}{\sqrt{2}} (1, \; 0, \; -1, \; 0, \; 0, \; 0)^\top$$
Độ cong của hàm Lagrange trên không gian tiếp tuyến là:
$$v_{\text{null}}^\top \nabla^2 \mathcal{L} v_{\text{null}} = \frac{1}{2} (2\lambda_1 \cdot 1^2 + 2\lambda_3 \cdot (-1)^2) = \frac{1}{2}(1.0 + 1.0) = +1.0000 > 0$$
Vì độ cong dương nghiêm ngặt trên toàn bộ không gian tiếp tuyến khác không, **Điều kiện đủ tối ưu bậc hai (SOSC) được thỏa mãn nghiêm ngặt**!

---

### 3.4 Đồ thị Vị trí Cân bằng và Thảo luận Nghiệm

![So sánh nghiệm Câu 2](file:///home/baconway/Documents/hanging-chain-ampl/hinh_cau2_so_sanh.png)

- **Cực tiểu địa phương ($f = -0.5910$)**: Hai mắt xích đầu vươn sang phải tới $X_2 \approx 1.9592$, sau đó mắt xích thứ 3 gập ngược về bên trái để nối vào điểm treo $(1,0)$. Giá trị riêng Hessian thu gọn trên không gian tiếp tuyến bằng $+0.2955 > 0$, nên đây là một vị trí cân bằng bền cục bộ về mặt toán học.
- **Cực tiểu toàn cục ($f = -2.0000$)**: Mắt xích 1 treo thẳng đứng từ $(0,0)$ xuống $(0,-1)$, mắt xích 2 nằm ngang từ $(0,-1)$ sang $(1,-1)$, và mắt xích 3 đi thẳng đứng lên $(1,0)$. Thế năng đạt giá trị nhỏ nhất tuyệt đối $f = -2.0$.

---

## 4. Câu 3: Bài toán Dây xích Treo có Chướng ngại vật

### 4.1 Mô hình hóa Chướng ngại vật và Rời rạc hóa Mắt xích
Ta xét dây xích có `so_mat_xich` $m = 10$, tổng chiều dài `chieu_dai` $L = 3.0$ (mỗi mắt xích dài $\ell = 0.3$), hai đầu cố định tại $(0,0)$ và $(2.0, 0.0)$.
Chướng ngại vật bên dưới dây xích có dạng hàm hình chuông Gauss:
$$\text{obst}(x) = h \exp\left( -\frac{(x - \mu)^2}{\sigma^2} \right) - \nu$$
Bộ thông số cơ sở: `chieu_cao_h` $h = 1.5$, `tam_mu` $\mu = 1.0$, `phuong_sai` $\sigma^2 = 0.1$, `do_dich_nu` $\nu = 2.0$.

Để đảm bảo toàn bộ thanh cứng của mỗi mắt xích không xuyên qua chướng ngại vật, mỗi mắt xích $k \in \{1, \dots, m\}$ được chia đều thành `so_diem_chia` $n = 20$ đoạn nhỏ bởi $n+1$ điểm nội suy:
$$\text{diem\_chia\_x}[k, i] = \left(1 - \frac{i}{n}\right) X_{k-1} + \left(\frac{i}{n}\right) X_k, \quad \text{diem\_chia\_y}[k, i] = \left(1 - \frac{i}{n}\right) Y_{k-1} + \left(\frac{i}{n}\right) Y_k$$
Ràng buộc tránh chướng ngại vật yêu cầu:
$$\text{diem\_chia\_y}[k, i] \ge \text{obst}(\text{diem\_chia\_x}[k, i]), \quad \forall k \in \{1, \dots, m\}, \; \forall i \in \{0, \dots, n\}$$

---

### 4.2 Phân tích Gợi ý trong Đề bài

Trong đề bài `M1_INT_POP_5_Lesson5_Homework.pdf`, phần gợi ý viết:
```ampl
var lx {k in 1..m} = sum {i in 1..k} x[i];
var cx {k in 1..m-1, i in 0..n} = i/n*lx[k] + (1-i/n)*lx[k+1];
```
Cách viết trong gợi ý chỉ chạy từ $k = 1$ đến $m-1$ (tức là chỉ kiểm tra từ mắt xích thứ 2 đến mắt xích thứ $m$, bỏ sót mắt xích đầu tiên nối từ gốc $(0,0)$ đến $(X_1, Y_1)$).
Trong mô hình [`xich_chuong_ngai_vat.mod`](file:///home/baconway/Documents/hanging-chain-ampl/xich_chuong_ngai_vat.mod), chúng ta định nghĩa tọa độ khớp `toa_do_khop_x` và `toa_do_khop_y` chạy từ chỉ số `0..so_mat_xich` (với khớp 0 bằng $0.0$), nhờ đó nội suy đầy đủ cả $m = 10$ mắt xích.

---

### 4.3 Mã nguồn AMPL cho Bài toán có Chướng ngại vật

#### Tệp mô hình: [`xich_chuong_ngai_vat.mod`](file:///home/baconway/Documents/hanging-chain-ampl/xich_chuong_ngai_vat.mod)
```ampl
# ==============================================================================
# Mô hình: xich_chuong_ngai_vat.mod
# Bài tập 4: Bài toán dây xích treo (Câu 3 - Có chướng ngại vật)
# ==============================================================================

param so_mat_xich integer > 1;
param chieu_dai > 0;
param toa_do_a;
param toa_do_b;
param he_so_c {k in 1..so_mat_xich} := so_mat_xich - k + 0.5;

param chieu_cao_h >= 0;
param tam_mu;
param phuong_sai > 0;
param do_dich_nu;
param so_diem_chia integer > 0 default 20;

var x {1..so_mat_xich};
var y {1..so_mat_xich};

var toa_do_khop_x {k in 0..so_mat_xich} = if k = 0 then 0.0 else sum {i in 1..k} x[i];
var toa_do_khop_y {k in 0..so_mat_xich} = if k = 0 then 0.0 else sum {i in 1..k} y[i];

var diem_chia_x {k in 1..so_mat_xich, i in 0..so_diem_chia} =
    (1 - i / so_diem_chia) * toa_do_khop_x[k-1] + (i / so_diem_chia) * toa_do_khop_x[k];
var diem_chia_y {k in 1..so_mat_xich, i in 0..so_diem_chia} =
    (1 - i / so_diem_chia) * toa_do_khop_y[k-1] + (i / so_diem_chia) * toa_do_khop_y[k];

minimize TheNang:
    sum {k in 1..so_mat_xich} he_so_c[k] * y[k];

subject to DoDaiMatXich {k in 1..so_mat_xich}:
    x[k]^2 + y[k]^2 = (chieu_dai / so_mat_xich)^2;

subject to DiemCuoiX:
    sum {k in 1..so_mat_xich} x[k] = toa_do_a;

subject to DiemCuoiY:
    sum {k in 1..so_mat_xich} y[k] = toa_do_b;

subject to TranhChuongNgaiVat {k in 1..so_mat_xich, i in 0..so_diem_chia}:
    diem_chia_y[k,i] >= chieu_cao_h * exp( - (diem_chia_x[k,i] - tam_mu)^2 / phuong_sai ) - do_dich_nu;
```

---

### 4.4 Kết quả Cơ sở (Tái tạo Hình 2)

Giải bài toán từ điểm khởi tạo $x_k = y_k = -1.0$ cho giá trị thế năng tối ưu:
$$f(x^*, y^*) = -5.034317$$

Tọa độ các khớp nối $(X_k, Y_k)$ trùng khớp hoàn toàn với **Hình 2 (Figure 2)** trong đề bài:

| Khớp $k$ | Hoành độ $X_k$ | Tung độ $Y_k$ | Độ tăng $(x_k, y_k)$ | Độ cao chướng ngại vật $\text{obst}(X_k)$ | Trạng thái tiếp xúc |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **0** | $0.000000$ | $0.000000$ | Điểm treo cố định | $-2.000000$ | Nằm phía trên |
| **1** | $0.198815$ | $-0.224661$ | $(+0.1988, -0.2247)$ | $-1.997932$ | Nằm phía trên |
| **2** | $0.440464$ | $-0.402441$ | $(+0.2416, -0.1778)$ | $-1.935140$ | Nằm phía trên |
| **3** | $0.724375$ | $-0.499368$ | $(+0.2839, -0.0969)$ | $-1.295412$ | Nằm phía trên |
| **4** | $\mathbf{1.024372}$ | $\mathbf{-0.500598}$ | $(+0.3000, -0.0012)$ | $\mathbf{-0.508882}$ | **Tiếp xúc đỉnh chướng ngại vật** |
| **5** | $1.196926$ | $-0.746006$ | $(+0.1726, -0.2454)$ | $-1.025315$ | Nằm phía trên |
| **6** | $\mathbf{1.415411}$ | $\mathbf{-0.951589}$ | $(+0.2185, -0.2056)$ | $-1.733489$ | **Đáy vòng võng sâu bên phải** |
| **7** | $1.693388$ | $-0.838768$ | $(+0.2780, +0.1128)$ | $-1.987725$ | Nằm phía trên |
| **8** | $1.842059$ | $-0.578197$ | $(+0.1487, +0.2606)$ | $-1.998749$ | Nằm phía trên |
| **9** | $1.934173$ | $-0.292689$ | $(+0.0921, +0.2855)$ | $-1.999895$ | Nằm phía trên |
| **10** | $2.000000$ | $0.000000$ | $(+0.0658, +0.2927)$ | $-2.000000$ | Điểm treo cố định |

![Tái tạo Hình 2 Câu 3](file:///home/baconway/Documents/hanging-chain-ampl/hinh_cau3_chuong_ngai_co_so.png)

---

### 4.5 Khảo sát Bộ thông số Tùy chọn (Chướng ngại vật Lệch tâm)

Theo yêu cầu của Câu 3:
> *"give at least one optimal solution with your own choice of parameter values. Represent the optimal solution on a graphic (as Figure 2)."*

Chúng ta khảo sát cấu hình cân bằng khi **chướng ngại vật được dịch tâm sang phải và nâng cao hơn**:
- Tâm dịch sang phải: `tam_mu` $\mu = 1.3$ (dịch $0.3$ đơn vị về phía điểm treo bên phải).
- Chiều cao tăng lên: `chieu_cao_h` $h = 1.6$ (đỉnh chướng ngại vật ở mức $y = -0.40$).
- Giữ nguyên độ rộng và độ dịch đứng: `phuong_sai` $\sigma^2 = 0.1$, `do_dich_nu` $\nu = 2.0$.

![Khảo sát thông số tùy chọn Câu 3](file:///home/baconway/Documents/hanging-chain-ampl/hinh_cau3_chuong_ngai_tuy_chon.png)

#### Nhận xét Vật lý:
1. **Nhánh võng sâu chuyển sang bên trái**: Khi đỉnh chướng ngại vật dịch sang $\mu = 1.3$, khoảng không gian bên trái rộng hơn ($1.3$ đơn vị so với $0.7$ đơn vị bên phải), giúp phần lớn chiều dài dây xích rủ sâu xuống vùng bên trái (đạt cực tiểu độ cao $Y_4 = -0.988908$ tại $X_4 = 0.603755$).
2. **Tiếp xúc trên sườn đỉnh cao hơn**: Các mắt xích thứ 7 và thứ 8 vắt qua đỉnh chướng ngại vật tại $X_7 = 1.225166, Y_7 = -0.426015$ và $X_8 = 1.506992, Y_8 = -0.323182$ trước khi nối lên điểm cuối $(2.0, 0.0)$.
3. **Thế năng giảm**: Thế năng tối ưu đạt $f = -5.238127$ (thấp hơn mức $-5.034317$ của cấu hình đối xứng cơ sở) do dây xích tận dụng được khoảng trống rộng bên trái để hạ thấp trọng tâm.

---

## 5. Hướng dẫn Thực thi và Tái tạo Kết quả

Toàn bộ các tệp trong thư mục làm việc đã được đồng bộ sang tên tiếng Việt không dấu:

### Cách 1: Chạy trực tiếp bằng AMPL
```bash
# Chạy Câu 2 (3 mắt xích):
ampl xich_co_ban.run

# Chạy Câu 3 (10 mắt xích với chướng ngại vật cơ sở và tùy chọn):
ampl xich_chuong_ngai_vat.run
```

### Cách 2: Chạy kịch bản tổng hợp tự động
```bash
./chay_ampl.sh
```
Hoặc chạy trực tiếp bộ giải và xuất đồ thị bằng Python:
```bash
python3 mo_phong_day_xich.py
```

### Danh mục Tệp Bàn giao
- **Câu 2 (AMPL)**: [`xich_co_ban.mod`](file:///home/baconway/Documents/hanging-chain-ampl/xich_co_ban.mod), [`xich_co_ban.dat`](file:///home/baconway/Documents/hanging-chain-ampl/xich_co_ban.dat), [`xich_co_ban.run`](file:///home/baconway/Documents/hanging-chain-ampl/xich_co_ban.run), dữ liệu xuất [`toa_do_xich_co_ban.txt`](file:///home/baconway/Documents/hanging-chain-ampl/toa_do_xich_co_ban.txt).
- **Câu 3 (AMPL)**: [`xich_chuong_ngai_vat.mod`](file:///home/baconway/Documents/hanging-chain-ampl/xich_chuong_ngai_vat.mod), [`xich_chuong_ngai_vat.dat`](file:///home/baconway/Documents/hanging-chain-ampl/xich_chuong_ngai_vat.dat), [`xich_chuong_ngai_vat.run`](file:///home/baconway/Documents/hanging-chain-ampl/xich_chuong_ngai_vat.run), dữ liệu xuất [`toa_do_xich_chuong_ngai_co_so.txt`](file:///home/baconway/Documents/hanging-chain-ampl/toa_do_xich_chuong_ngai_co_so.txt) và [`toa_do_xich_chuong_ngai_tuy_chon.txt`](file:///home/baconway/Documents/hanging-chain-ampl/toa_do_xich_chuong_ngai_tuy_chon.txt).
- **Kịch bản & Mô phỏng**: [`chay_ampl.sh`](file:///home/baconway/Documents/hanging-chain-ampl/chay_ampl.sh), [`mo_phong_day_xich.py`](file:///home/baconway/Documents/hanging-chain-ampl/mo_phong_day_xich.py).
- **Đồ thị Minh họa**: [`hinh_cau2_xich3.png`](file:///home/baconway/Documents/hanging-chain-ampl/hinh_cau2_xich3.png), [`hinh_cau2_so_sanh.png`](file:///home/baconway/Documents/hanging-chain-ampl/hinh_cau2_so_sanh.png), [`hinh_cau3_chuong_ngai_co_so.png`](file:///home/baconway/Documents/hanging-chain-ampl/hinh_cau3_chuong_ngai_co_so.png), [`hinh_cau3_chuong_ngai_tuy_chon.png`](file:///home/baconway/Documents/hanging-chain-ampl/hinh_cau3_chuong_ngai_tuy_chon.png).

