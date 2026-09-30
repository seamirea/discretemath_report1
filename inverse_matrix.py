from fractions import Fraction


class SingularMatrixError(ValueError):
    pass


def print_matrix(matrix, title=None):
    if title:
        print(title)
    widths = []
    for j in range(len(matrix[0])):
        width = 0
        for i in range(len(matrix)):
            width = max(width, len(str(matrix[i][j])))
        widths.append(width)
    for row in matrix:
        print("  [ ", end="")
        for j in range(len(row)):
            if j > 0:
                print("  ", end="")
            print(str(row[j]).rjust(widths[j]), end="")
        print(" ]")
    print()


def read_matrix():
    while True:
        try:
            n = int(input("Enter matrix size n: "))
            if n <= 0:
                raise ValueError
            break
        except ValueError:
            print("Error: n must be a positive integer.")

    print("Enter matrix:")
    matrix = []
    for i in range(n):
        while True:
            try:
                values = input(f"Row {i + 1}: ").split()
                if len(values) != n:
                    print(f"Error: enter exactly {n} numbers.")
                    continue
                row = []
                for value in values:
                    row.append(Fraction(value))
                matrix.append(row)
                break
            except (ValueError, ZeroDivisionError):
                print("Error: enter valid numbers. A fraction cannot have a zero denominator.")
    return matrix


def minor(matrix, excluded_row, excluded_column):
    result = []
    for i in range(len(matrix)):
        if i == excluded_row:
            continue
        row = []
        for j in range(len(matrix[i])):
            if j != excluded_column:
                row.append(matrix[i][j])
        result.append(row)
    return result


def determinant(matrix):
    n = len(matrix)
    if n == 0:
        return Fraction(1)
    if n == 1:
        return Fraction(matrix[0][0])
    if n == 2:
        return (Fraction(matrix[0][0]) * matrix[1][1]
                - Fraction(matrix[0][1]) * matrix[1][0])
    det = Fraction(0)
    for j in range(n):
        smaller_matrix = minor(matrix, 0, j)
        det += (-1) ** j * Fraction(matrix[0][j]) * determinant(smaller_matrix)
    return det


def inverse_by_determinant(matrix, show_steps=False):
    det = determinant(matrix)
    if show_steps:
        print(f"Step 1. det(A) = {det}\n")
    if det == 0:
        raise SingularMatrixError("det(A) = 0. No inverse exists.")

    n = len(matrix)
    cofactors = []
    for i in range(n):
        row = []
        for j in range(n):
            smaller_matrix = minor(matrix, i, j)
            cofactor = (-1) ** (i + j) * determinant(smaller_matrix)
            row.append(cofactor)
        cofactors.append(row)

    adjugate = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(cofactors[j][i])
        adjugate.append(row)

    result = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(adjugate[i][j] / det)
        result.append(row)
    if show_steps:
        print_matrix(cofactors, "Step 2. Cofactor matrix:")
        print_matrix(adjugate, "Step 3. Adjugate matrix (transpose of cofactor matrix):")
        print_matrix(result, "Step 4. A^(-1) = adj(A) / det(A):")
    return result


def identity_matrix(n):
    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append(Fraction(1))
            else:
                row.append(Fraction(0))
        matrix.append(row)
    return matrix


def inverse_by_gauss_jordan(matrix, show_steps=False):
    n = len(matrix)
    identity = identity_matrix(n)
    augmented = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(Fraction(matrix[i][j]))
        for j in range(n):
            row.append(identity[i][j])
        augmented.append(row)
    if show_steps:
        print_matrix(augmented, "Start with [A | I]:")

    for column in range(n):
        pivot_row = column
        for row in range(column + 1, n):
            if abs(augmented[row][column]) > abs(augmented[pivot_row][column]):
                pivot_row = row
        if augmented[pivot_row][column] == 0:
            raise SingularMatrixError(
                f"No nonzero pivot in column {column + 1}. "
                "No inverse exists.")
        if pivot_row != column:
            augmented[column], augmented[pivot_row] = augmented[pivot_row], augmented[column]
            if show_steps:
                print_matrix(augmented, f"Swap R{column + 1} and R{pivot_row + 1}:")

        pivot = augmented[column][column]
        for j in range(2 * n):
            augmented[column][j] = augmented[column][j] / pivot
        if show_steps:
            print_matrix(augmented, f"R{column + 1} <- R{column + 1} / ({pivot}):")

        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor == 0:
                continue
            for j in range(2 * n):
                augmented[row][j] -= factor * augmented[column][j]
            if show_steps:
                print_matrix(augmented,
                             f"R{row + 1} <- R{row + 1} - ({factor}) * R{column + 1}:")
    if show_steps:
        print_matrix(augmented, "Final augmented matrix [I | A^(-1)]:")
    result = []
    for i in range(n):
        row = []
        for j in range(n, 2 * n):
            row.append(augmented[i][j])
        result.append(row)
    return result


def multiply_matrices(left, right):
    result = []
    for i in range(len(left)):
        row = []
        for j in range(len(right[0])):
            total = Fraction(0)
            for k in range(len(right)):
                total += left[i][k] * right[k][j]
            row.append(total)
        result.append(row)
    return result


def main():
    print("Inverse Matrix Calculator")
    print("Method 1: Determinant / Adjugate")
    print("Method 2: Gauss-Jordan Elimination\n")
    print("Use small matrices for easy-to-read calculation steps.")
    print("The determinant method can be slow for large matrices.\n")
    matrix = read_matrix()
    print_matrix(matrix, "\nInput matrix A:")
    results = []
    methods = [("Method 1: Determinant / Adjugate", inverse_by_determinant),
               ("Method 2: Gauss-Jordan Elimination", inverse_by_gauss_jordan)]
    for title, method in methods:
        print(title)
        try:
            result = method(matrix, show_steps=True)
            print_matrix(result, "Inverse matrix:")
            results.append(result)
        except SingularMatrixError as error:
            print(f"Error: {error}\n")
            results.append(None)

    print("\n=== Comparison of Two Methods ===")
    first, second = results
    if first is None and second is None:
        print("Both methods found that no inverse exists.")
    elif first is None or second is None:
        print("Results do not match: only one method found an inverse.")
    else:
        if first == second:
            print("Results match.")
        else:
            print("Results do not match.")
        print("\n=== Additional Feature: Inverse Verification ===")
        print("Check whether A * A^(-1) = I\n")
        for label, inverse in [("Determinant / Adjugate", first), ("Gauss-Jordan", second)]:
            product = multiply_matrices(matrix, inverse)
            print_matrix(product, f"{label}: A * A^(-1):")
            if product == identity_matrix(len(matrix)):
                print("Check passed: the result is the identity matrix.\n")
            else:
                print("Check failed.\n")


if __name__ == "__main__":
    try:
        main()
    except (EOFError, KeyboardInterrupt):
        print("\nInput stopped.")
