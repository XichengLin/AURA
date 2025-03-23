function angles = compute_angles(P1, P2)
    % 计算 P2 相对于 P1 的方位角和俯仰角
    % P1 和 P2 均为 3×1 的列向量 [x; y; z]
    % 返回 angles 也是 2×1 列向量 [azimuth; elevation]

    % 计算坐标差值
    delta = P2 - P1;

    % 提取坐标差分量
    delta_x = delta(1);
    delta_y = delta(2);
    delta_z = delta(3);

    % 计算方位角（azimuth）
    azimuth = atan2d(delta_y, delta_x);

    % 计算俯仰角（elevation），单位为度
    elevation = atan2d(delta_z, sqrt(delta_x^2 + delta_y^2));

    % 以列向量返回
    angles = [azimuth; elevation];
end