import { defineConfig } from 'vitest/config'

export default defineConfig({
  test: {
    include: ['release-check/**/*.test.ts'],
    testTimeout: 3_600_000,
  },
})
