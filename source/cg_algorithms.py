import math


def draw_line(p_list, algorithm):
    """绘制线段

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """

    "对空集的处理"
    if len(p_list) == 0:
        return []

    x0, y0 = p_list[0]
    x1, y1 = p_list[1]

    "对相同的起始点和终末点的处理"
    if [x0, y0] == [x1, y1]:
        return [[x0, y0]]

    result = []
    if algorithm == 'DDA':
        dx = x1 - x0
        dy = y1 - y0
        if(abs(dx) > abs(dy)):
            steps = abs(dx)
        else:
            steps = abs(dy)
        delta_x = dx / steps
        delta_y = dy / steps
        for i in range(steps+1):
            result.append([int(x0+i*delta_x), int(y0+i*delta_y)])
        pass
    elif algorithm == 'Bresenham':
        dx = abs(x1 - x0)
        dy = abs(y1 - y0)
        sx1 = sx2 = 1 if x0 < x1 else -1
        sy1 = sy2 = 1 if y0 < y1 else -1
        sx2 = sx1
        sy2 = sy1
        if(dx > dy):
            faststeps = dx
            slowsteps = dy
            sy2 = 0
        else:
            faststeps = dy
            slowsteps = dx
            sx2 = 0
        err = faststeps//2
        while True:
            result.append([x0, y0])
            if(x0 == x1 and y0 == y1):
                break
            err += slowsteps
            if(err > faststeps):
                err -= faststeps
                x0 += sx1
                y0 += sy1
            else:
                x0 += sx2
                y0 += sy2
        pass
    return result


def draw_polygon(p_list, algorithm):
    """绘制多边形

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 多边形的顶点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    for i in range(len(p_list)):
        line = draw_line([p_list[i - 1], p_list[i]], algorithm)
        result += line
    return result


def draw_ellipse(p_list):
    """绘制椭圆（采用中点圆生成算法）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 椭圆的矩形包围框左上角和右下角顶点坐标
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    if x0 > x1:
        x0, x1 = x1, x0
    if y0 > y1:
        y0, y1 = y1, y0
    result = []
    a = (x1 - x0)/2
    b = (y1 - y0)/2
    mx = int(x0 + a)
    my = int(y0 + b)
    sqa = a*a
    sqb = b*b
    x = int(0)
    y = int(b)
    p1 = sqb - sqa*b + 0.25*sqa
    while sqb*x < sqa*y:
        result.append([mx + x, my + y])
        result.append([mx + x, my - y])
        result.append([mx - x, my + y])
        result.append([mx - x, my - y])
        if p1 < 0:
            x = x+1
            p1 = p1 + 2*sqb*x + sqb
        else:
            x = x + 1
            y = y - 1
            p1 = p1 + 2*sqb*x-2*sqa*y+sqb
    p2 = sqb*(x+0.5)*(x+0.5)+sqa*(y-1)*(y-1)-sqa*sqb
    while y >= 0:
        result.append([mx + x, my + y])
        result.append([mx + x, my - y])
        result.append([mx - x, my + y])
        result.append([mx - x, my - y])
        if p2 > 0:
            y = y-1
            p2 = p2 - 2*sqa*y + sqa
        else:
            x = x + 1
            y = y-1
            p2 = p2 + 2*sqb*x - 2*sqa*y + sqa
    return result


def draw_curve_Bezier(p_list, n):
    def one_Bezier(x0, x1, t):
        return (1 - t) * x0 + t*x1

    def n_Bezier(xs, n, k, t):
        if n == 1:
            return one_Bezier(xs[k], xs[k + 1], t)
        else:
            return (1 - t) * n_Bezier(xs, n-1, k, t) + t*n_Bezier(xs, n-1, k+1, t)

    xs = []
    ys = []
    for i in range(n):
        xs += [p_list[i][0]]
        ys += [p_list[i][1]]
    result_xs = []
    result_ys = []
    t = 0.0
    step = 1 / (n * 5000)
    while t < 1.0:
        result_xs += [n_Bezier(xs, n - 1, 0, t)]
        result_ys += [n_Bezier(ys, n - 1, 0, t)]
        t += step
    result = []
    for i in range(len(result_xs)):
        result += [(int(result_xs[i]), int(result_ys[i]))]
    result = list(set(result))
    return result


def draw_curve_Bspline(p_list, n):
    def Bspline(xs, k, t):
        b0 = 1/6 * ((1-t)**3)
        b1 = 1/6 * (3*t**3 - 6*t**2 + 4)
        b2 = 1/6 * (-3*t**3 + 3*t**2 + 3*t + 1)
        b3 = 1/6 * (t**3)
        return b0*xs[k] + b1*xs[k+1] + b2 * xs[k+2] + b3 * xs[k+3]
    result = []

    xs = []
    ys = []
    for i in range(n):
        xs += [p_list[i][0]]
        ys += [p_list[i][1]]
    result_xs = []
    result_ys = []
    for k in range(0, n-3):
        t = 0.0
        step = 1 / 5000
        while t < 1.0:
            result_xs += [Bspline(xs, k, t)]
            result_ys += [Bspline(ys, k, t)]
            t += step
    result = []
    for i in range(len(result_xs)):
        result += [(int(result_xs[i]), int(result_ys[i]))]
    result = list(set(result))
    return result


def draw_curve(p_list, algorithm):
    """绘制曲线

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 曲线的控制点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'Bezier'和'B-spline'（三次均匀B样条曲线，曲线不必经过首末控制点）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    n = len(p_list)
    if algorithm == 'Bezier':
        result = draw_curve_Bezier(p_list, n)
        pass
    elif algorithm == 'B-spline':
        result = draw_curve_Bspline(p_list, n)
    print(result)
    return result


def translate(p_list, dx, dy):
    """平移变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param dx: (int) 水平方向平移量
    :param dy: (int) 垂直方向平移量
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    result = []
    for x, y in p_list:
        x += dx
        y += dy
        result.append([x, y])
    return result


def rotate(p_list, x, y, r):
    """旋转变换（除椭圆外）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 旋转中心x坐标
    :param y: (int) 旋转中心y坐标
    :param r: (int) 顺时针旋转角度（°）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    result = []
    radr = math.radians(r)
    for xi, yi in p_list:
        d0 = math.sqrt((xi-x)**2+(yi-y)**2)
        r0 = math.atan2(yi-y, xi-x)
        result.append([int(x + d0*math.cos(r0+radr)),
                       int(y + d0 * math.sin(r0+radr))])
    return result


def scale(p_list, x, y, s):
    """缩放变换
    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 缩放中心x坐标
    :param y: (int) 缩放中心y坐标
    :param s: (float) 缩放倍数
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    result = []
    for xi, yi in p_list:
        dx = xi - x
        dy = yi - y
        result.append([int(x + s*dx), int(y + s*dy)])
    return result


def clip_CohenSutherland(p_list, x_min, y_min, x_max, y_max):
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]

    LEFT, RIGHT, TOP, BOTTOM = 1, 2, 4, 8

    def encode(x, y):
        c = 0
        if x < x_min:
            c = c | LEFT
        elif x > x_max:
            c = c | RIGHT
        if y < y_min:
            c = c | TOP
        elif y > y_max:
            c = c | BOTTOM
        return c

    code0 = encode(x0, y0)
    code1 = encode(x1, y1)

    while True:
        if code0 == 0 and code1 == 0:
            "完全在视窗内：返回裁剪好的线段"
            return [[x0, y0], [x1, y1]]
        if code0 & code1 != 0:
            "完全不在视窗内：返回空集"
            return []

        "点[x0,y0]在视窗外，否则交换两点"
        if code0 == 0:
            x0, y0, x1, y1 = x1, y1, x0, y0
            code0, code1 = code1, code0

        "计算新点"
        if LEFT & code0 != 0:
            if x0 == x1:
                y0 = y1
            else:
                y0 = int(y0+(y1-y0)*(x_min - x0)/(x1-x0))
            x0 = x_min
        elif RIGHT & code0 != 0:
            if x0 == x1:
                y0 = y1
            else:
                y0 = int(y0+(y1-y0)*(x_max - x0)/(x1-x0))
            x0 = x_max
        elif TOP & code0 != 0:
            if y0 == y1:
                x0 = x1
            else:
                x0 = int(x0+(x1-x0)*(y_min - y0)/(y1-y0))
            y0 = y_min
        elif BOTTOM & code0 != 0:
            if y0 == y1:
                x0 = x1
            else:
                x0 = int(x0+(x1-x0)*(y_max - y0)/(y1-y0))
            y0 = y_max
        "对新点进行编码"
        code0 = encode(x0, y0)

        "while循环的结束"
        print([x0, y0], [x1, y1])
    return ([x0, y0], [x1, y1])


def clip_LiangBarsky(p_list, x_min, y_min, x_max, y_max):
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]

    def _clip(p, q, u) -> int:
        if(p < 0):
            r = q/p
            if(r > u[1]):
                return False
            if(r > u[0]):
                u[0] = r
                return True
        elif(p > 0):
            r = q / p
            if(r < u[0]):
                return False
            if(r < u[1]):
                u[1] = r
                return True
        else:
            return q >= 0

    u = [0, 1]
    dx = x1 - x0
    dy = y1 - y0

    if(_clip(-dx, x0-x_min, u)):
        if(_clip(dx, x_max - x0, u)):
            if(_clip(-dy, y0 - y_max, u)):
                if(_clip(dy, y_min - y0, u)):
                    nx0, ny0 = x0 + u[0]*dx, y0+u[0]*dy
                    nx1, ny1 = x0 + u[1]*dx, y0 + u[1]*dy
                    return ([nx0, ny0], [nx1, ny1])
        pass
    return []


def clip(p_list, x_min, y_min, x_max, y_max, algorithm):
    """线段裁剪
    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param x_min: 裁剪窗口左上角x坐标
    :param y_min: 裁剪窗口左上角y坐标
    :param x_max: 裁剪窗口右下角x坐标
    :param y_max: 裁剪窗口右下角y坐标
    :param algorithm: (string) 使用的裁剪算法，包括'Cohen-Sutherland'和'Liang-Barsky'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1]]) 裁剪后线段的起点和终点坐标
    """

    if algorithm == 'Cohen-Sutherland':
        return clip_CohenSutherland(p_list, x_min, y_min, x_max, y_max)
    elif algorithm == 'Liang-Barsky':
        return clip_LiangBarsky(p_list, x_min, y_min, x_max, y_max)
    return []
